# Libraries ---------------------------------------------------------------
library(readr)
library(magrittr)
library(dplyr)
library(chemhelper) # https://github.com/stanstrup/chemhelper


peaklist_mzmine_3 <- read_csv("Z:/JPZS - Jan Stanstrup/_Auto_annotation/mzMine_example.csv") %>% 
    select(-ncol(.)) %>% # last column seems empty. always?
    as.data.frame
# Read the peaklist into R ------------------------------------------------
peaklist_mzmine <- read_csv("Z:/_Data/LIP1/0034_Epos/Processed_Aug_2016/Overfeeding_Lipidomics_reannotated.csv") %>% 
    as.data.frame


# # Load lipidmaps db -----------------------------------------------------
lipidmaps <- readRDS("Z:/JPZS - Jan Stanstrup/lipid_maps_annotation/SDF_table.rds") %>% 
    as.data.frame(stringsAsFactors = FALSE) %>% 
    mutate(mz = as.numeric(EXACT_MASS)+1.0073) %>% 
    filter(!is.na(mz))



# Lets make adducts for all the compounds ---------------------------------
temp1 <- lipidmaps %>% mutate(mz = mz + 17.02655) %>% mutate(COMMON_NAME = paste0(COMMON_NAME," [+NH4+]"))
temp2 <- lipidmaps %>% mutate(mz = mz + 21.98194) %>% mutate(COMMON_NAME = paste0(COMMON_NAME," [+Na+]"))

lipidmaps %<>% bind_rows(temp1,temp2)
rm(temp1,temp2)



# Annotate dataset --------------------------------------------------------
assigned <- db.comp.assign(mz = peaklist_mzmine[,"row m/z"],
                           rt = peaklist_mzmine[,"row retention time"],
                           comp_name_db = lipidmaps$COMMON_NAME,
                           mz_db        = lipidmaps$mz,
                           rt_db        = rep(5,length(lipidmaps$mz)), # we just fake a list of retention times
                           mzabs        = 0.01,
                           ppm          = 10,
                           ret_tol      = Inf
) %>% 
    data.frame(anno_lipidmaps = ., stringsAsFactors = FALSE)


peaklist_mzmine %<>% bind_cols(assigned)
#rm(assigned)

Annotated_Lipid<- data.frame(aassigned$anno_lipidmaps,aassigned_2$anno_lipidmaps,peaklist_mzmine)
write.csv(Annotated_Lipid, "Z:/_Data/LIP1/0034_Epos/Processed_Aug_2016/Overfeeding_Lipidomics_reannotated_2.csv", row.names = FALSE)
