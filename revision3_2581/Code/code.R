# -------------------------------------------------------------------------
# Supplemental R Code to reproduce the results of
# Szöcs, Schäfer. webchem: An R Package to Retrieve Chemical Information
#   from the Web
# ------------------------------------------------------------------------


# Install webchem from CRAN -----------------------------------------------
#install.packages("webchem")


# load packages -----------------------------------------------------------
#library("webchem")
library("ggplot2")
library("testthat")

# load datasets -----------------------------------------------------------
# load jagst dataset
data("jagst", package = "webchem")
# print first 6 lines of data
head(jagst)

# lc50 dataset
data("lc50", package = "webchem")
head(lc50)


# Use Case 1: Query identifiers -------------------------------------------
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

# query SMILES from PubChem
cids <- get_cid(etox_cas)
pc_data <- pc_prop(cids, properties = "CanonicalSMILES")
pc_smiles <- smiles(pc_data)

# ChemSpider needs an API key.
# This is a webchem specific token.
# Please use it only for reproduction / testing.
token <- "37bf5e57-9091-42f5-9274-650a64398aaf"
# query InChiKey from ChemSpider
#csids <- get_csid(pc_smiles, from = "smiles", apikey = apikey)
#cs_inchikey <- cs_convert(csids, from = "csid", to = "inchikey", apikey = apikey)

# combine into single data.frame
#res <- data.frame(name = subs, cas = etox_cas, smiles = pc_smiles,
#                  cid = pc_data$CID, inchikey = cs_inchikey,
#                  csid = cs_data$csid,
#                  stringsAsFactors = FALSE)
#head(res)

# Automatic tests to test if results are reproduced correctly
test_that("Use Case 1: Query identifiers",{

  expect_equal(names(etox_cas), c("8932", "8494", NA, "8397", "7240", "7331"))
  expect_equal(unname(etox_cas), c("105-67-9", "1570-64-5", NA, "1912-24-9",
                                   "71-43-2", "6190-65-4"))
  expect_is(cids, "list")
  expect_length(cids, 6)

  expect_is(pc_data, c("pc_prop", "data.frame"))
  expect_dim(pc_data, c(6,2))

  expect_is(pc_smiles, "character")
  expect_length(pc_smiles, 6)

})


# Use Case II: Toxicity of different pesticide groups ---------------------
# query the pesticide compendium using CAS-numbers
aw_data <- aw_query(lc50$cas[1:3], type = "cas")

# shows internal structure of the data-object
str(aw_data[[1]])
# extract chemical group from list
igroup <- sapply(aw_data, function(y) y$subactivity[1])
igroup[1:3]

aw_data <- aw_query(lc50$cas, type = "cas")
igroup <- sapply(aw_data, function(y) unname(y["subactivity"][1]))


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

# Automatic tests to test if results are reproduced correctly
test_that("Use Case II: Toxicity of different pesticide groups",{

  expect_equal(names(igroup[1:3]), c("50-29-3", "52-68-6", "55-38-9"))
  expect_equal(unname(igroup[1:3]), c("organichlorine insecticides",
                                      "phosphonate insecticides",
                                      "phenyl organothiophosphate
                                      insecticides"))

})


# Use Case III: Name 10 most toxic chemicals ------------------------------
lc50_3 <- lc50[ order(lc50$value), ][1:3, ]
# query ChEBI names from CAS
lite <- get_chebiid(lc50_3$cas, match = 'best')
lite <- do.call(rbind, lite) # bind to data.frmae
lite$cas <- rownames(lite)
# query ChEBI complete entity
comp <- chebi_comp_entity(lite$chebiid) # use ChEBI ids to query the service
comp <- do.call(rbind, lapply(comp, `[[`, "parents")) # extract parents list and bind them to a data.frame
comp$chebiid <- gsub('\\.[0-9]+', '', rownames(comp))
comp <- comp[ comp$type == "has role", ] # refine to chemical roles
role <- aggregate(chebiName ~ chebiid, data = comp, FUN = function(x) paste0(x, collapse = ', '))
role <- merge(lite[ names(lite) %in% c('cas', 'chebiid', 'chebiasciiname') ], role, by = 'chebiid')
setNames(role, c('chebiid', 'name', 'cas', 'roles'))

test_that("Use Case III: Name 10 most toxic chemicals",{

})

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

test_that("Use Case IV: Querying partitioning coefficients",{
  expect_is(cid, "integer")
  expect_is(pc_data, c("pc_prop", "data.frame"))
  expect_is(lc50$logp, "numeric")
  expect_length(lc50$logp[!is.na(lc50$logp)], 118)
  expect_length(lc50$logp[is.na(lc50$logp)], 6)
  expect_equal(round(coef(mod),2), c(2.86, -0.36))
  expect_equal(round(summary(mod$sigma), 2), 1.47)
})

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

test_that("Use Case V: Regulatory information",{
  expect_is(mac, "data.frame")
  expect_equal(dim(mac), c(6,3))
  expect_equal(names(mac), c("etoxid","query","mac"))
  expect_equal(mac$etoxid, c(8397, 7240, 8836, 7442, 7571, 8756))
  expect_equal(mac$query, c("Atrazin", "Benzol", "Irgarol", "Isoproturon",
                            "Simazin", "Terbutryn"))
  expet_equal(mac$mac, c(2.000, 50.000, 0.016, 1.000, 4.000, 0.034))
  expect_is(jagst_eqs, "data.frame")
  expect_equal(names(jagst_eqs), c("substance", "date", "value", "qual",
                                   "etoxid", "mac"))
  expect_equal(head(jagst_eqs$substance), rep("Atrazin", times = 6))
  expect_equal(head(jagst_eqs$date), as.Date(c("2013-09-10", "2013-10-08",
                                               "2013-03-26", "2013-04-23",
                                               "2013-06-18", "2013-07-16")))
  expect_equal(head(jagst_eqs$value, c(0.0068, 0.0072, 0.0040, 0.0048, 0.0048,
                                       0.0052)))
  expect_equal(head(jagst_eqs$qual), rep("=", times = 6))
  expect_equal(head(jagst_eqs$etoxid), rep(8397, times = 6))
  expect_equal(head(jagst_eqs$mac), rep(2, times = 6))
})

# Utilitymac functions -------------------------------------------------------
# simple formatting check
is.inchikey("BQJCRHHNABKAKU-KBQPJGBKS-AN")
# formatting check
is.cas("64-17-6")

# or using ChemSpider
is.inchikey("BQJCRHHNABKAKU-KBQPJGBKSA-5", type = "chemspider")

test_that("Utility functions",{
  expect_false(is.inchikey("BQJCRHHNABKAKU-KBQPJGBKS-AN"))
  expect_false(is.cas("64-17-6"))
  expect_false(is.inchikey("BQJCRHHNABKAKU-KBQPJGBKSA-5", type = "chemspider"))
})

