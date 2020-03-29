# -------------------------------------------------------------------------
# Supplemental R Code to reproduce the results of Szöcs, Stirling, Scharmüller, Schäfer. 
# webchem: An R Package to Retrieve Chemical Information from the Web
# ------------------------------------------------------------------------

# Install webchem from CRAN -----------------------------------------------
library(devtools)
install_github("ropensci/webchem")

# load packages -----------------------------------------------------------
library(webchem)
library(ggplot2)

# load datasets -----------------------------------------------------------
# load jagst dataset
data("jagst", package = "webchem")
# print first 6 lines of data
head(jagst)

# lc50 dataset
data("lc50", package = "webchem")
head(lc50)

# Use Case I: Query identifiers -------------------------------------------
# unique substance names
subs <- unique(jagst$substance)

# search ETOX IDs, keeping only the best match
ids <- get_etoxid(subs, match = "best")
head(ids)

# use this id to query information from ETOX
etox_data <- etox_basic(ids$etoxid)
# extract CAS numbers from retrieved data
etox_cas <- cas(etox_data)
head(etox_cas)

# query other identifiers from other resources
# get_cid() returned the PubChem ID of sodium when the query was NA. This was
# fixed in the dev version.
# query SMILES from PubChem
cids <- get_cid(etox_cas)
pc_data <- pc_prop(cids, properties = "CanonicalSMILES")
pc_smiles <- smiles(pc_data)

# ChemSpider needs an API key and RSC explicitly asks users to keep the key
# secure. In order to test these functions, please register for an API key.
csids <- get_csid(pc_smiles, from = "smiles", apikey = apikey)
cs_inchikey <- cs_convert(csids$csid, from = "csid", to = "inchikey",
                          apikey = apikey)

# combine into single data.frame
res <- data.frame(name = subs,
                 cas = etox_cas,
                 smiles = pc_smiles,
                 cid = pc_data$CID,
                 inchikey = cs_inchikey,
                 csid = csids$csid,
                 stringsAsFactors = FALSE,
                 row.names = NULL)
head(res)

# Use Case II: Toxicity of different pesticide groups ---------------------
# query the pesticide compendium using CAS-numbers
aw_data <- aw_query(lc50$cas, type = "cas")

# shows internal structure of the data-object
str(aw_data[[1]])
# extract chemical group from list
igroup <- sapply(aw_data, function(y) {
  if (!is.na(y)) {
    y$subactivity[1]  
  } else {
    NA
  }
})
igroup[1:3]

# cleanup
lc50$type <- ifelse(grepl("carbamat", igroup), "Carbamates",
              ifelse(grepl("neonicotinoid", igroup), "Neo-\nnicotinoids",
                ifelse(grepl("pyrethroid", igroup), "Pyrethroids",
                  ifelse(grepl("organophosphate", igroup), "Organo-\nphosphates",
                         "other"))))
lc50$type <- factor(lc50$type, levels = c("Pyrethroids", "Carbamates",
                                          "Organo-\nphosphates",
                                          "Neo-\nnicotinoids",
                                          "other"))

# plot
p <- ggplot(lc50, aes(x = type, y = value)) +
  geom_boxplot(fill = "grey75") +
  geom_jitter(width = 0.5) +
  scale_y_log10() +
  labs(y = expression(LC[50~"D.magna, 48h"]), x = "") +
  theme_bw() +
  theme(text = element_text(size = 16))
p

# Use Case III: Name 10 most toxic chemicals ------------------------------
cas_rns <- lc50[order(lc50$value)[1:3],"cas"]
chebiids <- get_chebiid(cas_rns)
comp <- chebi_comp_entity(chebiids$chebiid)
pars <- lapply(comp, function(x) with(x, parents[parents$type == "has role",]))

# Use Case IV: Querying partitioning coefficients -------------------------
# query PubChem DATABASE
cid <- get_cid(lc50$cas, first = TRUE)
pc_data <- pc_prop(cid)
# lookat internal structure
str(pc_data)
# extract logP from properties data.frame
lc50$logp <- pc_data$XLogP

# model
mod <- lm(log10(value) ~ logp, data = lc50)
coef(mod) # coeficients
summary(mod)$sigma # RMSE

# plot
p <- ggplot(lc50, aes(x = logp, y = value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  scale_x_continuous(breaks = c(0, 2, 4, 6),
                     labels = format(10^c(0, 2, 4, 6),
                                     scientific = FALSE)) +
  labs(x = "P", y = expression(LC[50])) +
  theme_bw()
p

# Use Case V: Regulatory information --------------------------------------
# search EQS from ETOX using the already queried ETOX_IDs
eqs <- etox_targets(ids$etoxid)
# extract only MAC-EQS for the EU from the results
ids$mac <- sapply(eqs, function(y){
  if (length(y) == 1 && is.na(y)) {
    return(NA)
  } else {
    res <- y$res
    min(res[res$Country_or_Region == "EEC / EU" &
              res$Designation == "MAC-EQS", "Value_Target_LR"])
  }
})
# keep only compounds with MAC-EQS
(mac <- with(ids, ids[!is.na(mac) & is.finite(mac),
                      c("etoxid", "query", "mac")]))

# join with original data
jagst_eqs <- merge(jagst, mac, by.x = "substance", by.y = "query")
head(jagst_eqs)

# Utility functions -------------------------------------------------------
# simple formatting check
is.inchikey("BQJCRHHNABKAKU-KBQPJGBKS-AN")
# formatting check
is.cas("64-17-6")

# or using ChemSpider
is.inchikey("BQJCRHHNABKAKU-KBQPJGBKSA-5", type = "chemspider")

