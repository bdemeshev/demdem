library("readr")

parallel <- TRUE

create_model_list <- function() {
  df <- expand.grid(type="mcmc", chain_no=1:5, status="not estimated", 
                  sample=c("west","east","all"), W=c("border","road"))
  df <- df %>% mutate_each("as.factor",type, status, sample, W) %>% mutate(seed=chain_no+42)
  df <- df %>% mutate(file=paste0(type,"_",W,"_",sample,"_",chain_no,".Rds"), id=row_number())
return(df)
}

# create model list
mlist <- create_model_list()
write_csv(mlist, path = "../estimation/mlist_A.csv")

# estimate models from list
mlist <- read_csv("../estimation/mlist_A.csv")
mlist <- estimate_models(mlist, parallel = parallel) # status and filename are updated
write_csv(mlist, path = "../estimation/mlist_A.csv")

# create prediction list
# plist <- create_prediction_list(...)
# write_csv(plist, path = "../estimation/plist_A.csv")

# make predictions
# plist <- read_csv("../estimation/plist_A.csv")
# plist <- make_predictions(plist, parallel = parallel) # status and filename are updated
# write_csv(plist, path = "../estimation/plist_A.csv")

# evaluate models
# plist <- read_csv("../estimation/plist_A.csv")
# evlist <- evaluate_models(plist)



