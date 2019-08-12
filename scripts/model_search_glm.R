rm(list = ls())

.packages = c("dplyr", "tibble", "boot", "Rmisc", "purrr", "reshape2", "parallel", "optparse", "pracma", "broom")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only=TRUE)


option_list = list(
  make_option(c("-d", "--data"), type="character", default=NULL,
              help="input data", metavar="character"),
  make_option(c("-v", "--variabley"), type="character", default=NULL,
              help="name of the dependent variable as given in the input data", metavar="character"),
  make_option(c("-p", "--npredictors"), type="integer", default=2,
              help="max number of predictors in a glm [default: %default]", metavar="integer"),
  make_option(c("-o", "--output"), type = "character", default = "model_search_2wayglm_out.csv",
              help = "name of the output file [default: %default]", metavar = "character"))

opt_parser = OptionParser(option_list=option_list,
                          description = "Perform a series of generalised linear models (family binomial) with 2-way interactions for all possible combinations of an established set of predictors.",
                          epilogue = "Example: Rscript scripts/model_search_glm.R -d data/CZ_all_mating_clean.csv -v mountYNcontact -p 5 -o tables/model_search/glm_2way_5predictors.csv")
opt = parse_args(opt_parser)

if (is.null(opt$data) | is.null(opt$variabley)){
  print_help(opt_parser)
  stop("At least two arguments must be supplied (input data and name of the dependent variable).\n", call.=FALSE)
}

# outfile = "tables/model_search/model_search_2wayglm_out.csv"
dir.create(dirname(opt$output))
CZ_data = read.csv(opt$data, sep = ";")
# CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
# colnames(CZ_data)
CZ.glm = CZ_data[, c("log_female", "shape", "size_ratio", "size_ratio2", "ref_ecotype", "shore", "mountYNcontact")]
CZ.glm$size_ratio3 = CZ.glm$size_ratio^3
CZ.glm$ref_ecotype=as.integer(CZ.glm$ref_ecotype) # 1 for crab and 2 for wave
CZ.glm$shore=as.integer(CZ.glm$shore) # 1 for CZA, 2 for CZB, 3 for CZC, 4 for CZD
CZ.glm$eco_fem_int=CZ.glm$ref_ecotype*CZ.glm$log_female
CZ.glm$eco_shape_int=CZ.glm$ref_ecotype*CZ.glm$shape
CZ.glm$eco_ratio_int=CZ.glm$ref_ecotype*CZ.glm$size_ratio
CZ.glm$eco_ratio2_int=CZ.glm$ref_ecotype*CZ.glm$size_ratio2
CZ.glm$eco_ratio3_int=CZ.glm$ref_ecotype*CZ.glm$size_ratio3
CZ.glm$eco_shore_int=CZ.glm$ref_ecotype*CZ.glm$shore
CZ.glm$shape_fem_int=CZ.glm$shape*CZ.glm$log_female
CZ.glm$shape_ratio_int=CZ.glm$shape*CZ.glm$size_ratio
CZ.glm$shape_ratio2_int=CZ.glm$shape*CZ.glm$size_ratio2
CZ.glm$shape_ratio3_int=CZ.glm$shape*CZ.glm$size_ratio3
CZ.glm$shape_shore_int=CZ.glm$shape*CZ.glm$shore
CZ.glm$shore_fem_int=CZ.glm$shore*CZ.glm$log_female
CZ.glm$shore_ratio_int=CZ.glm$shore*CZ.glm$size_ratio
CZ.glm$shore_ratio2_int=CZ.glm$shore*CZ.glm$size_ratio2
CZ.glm$shore_ratio3_int=CZ.glm$shore*CZ.glm$size_ratio3
# colnames(CZ.glm)
# yvar_id = which(names(CZ.glm)=="mountYNcontact")
yvar_id = which(names(CZ.glm)==opt$variabley)
cat("\nSaving additive effects ...\n")
CZ_form_add = names(CZ.glm)[-yvar_id]
cat("\nSaving 2-way interaction effects ...\n")
CZ_form_int = unlist(lapply(2, function(n) combn(CZ_form_add, n, FUN=function(row) paste0(row, collapse = ":"))))
X = c(CZ_form_add, CZ_form_int)
# X = X[1:6]
cat("\nCreating models with additive effects and 2-way interactions ...\n")
CZ_form = unlist(lapply(1:5, function(n) {
  combn(X, n, FUN=function(row) paste0(opt$variabley, " ~ ", paste0(row, collapse = "+")))
}))
cat("\nRunning", length(CZ_form), "GLMs family binomial ...\n")
CZ_crit = bind_rows(lapply(CZ_form, function(frml) {
  a = glance(glm(frml,family = binomial, data=CZ.glm))
  a$frml = frml
  return(a)
}))
cat("\nSaving summary stats in", opt$output, "...\n")
write.csv(x = CZ_crit, file = opt$output, row.names = FALSE)
