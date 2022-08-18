
library(tidyr)
library(dplyr)
library(tidyverse)
library(writexl)
library(ggfortify)
library(FactoMineR)
library(factoextra)
library(MASS)
library(reshape2)
library(vegan)
library(ggrepel)
library(janitor)
library(ssh)






###Automating the file Transfer
#################### THIS IS YOUR INPUT/OUTPUT SECTION #####################

biobake_input_folder <- "~/biobakery"

Metaphlan_analysis_output_folder <- "~/VIRGO_Valencia_Results_K_shotgun/VIRGO_Valencia_Input/"

Metaphlan_analysis_outputfilename <- "K_shotgun_metaphlan_fsingle.csv"

Valencia_output <- "K_shotgun_metaphlan_CST_fsingle.csv"

Valencia_output_file_destination <- "VIRGO_Valencia_Results_K_shotgun/Valencia"



xlabels2=c("K12D04" = "-11", "K12D06" = "-9", "K12D09" = "-6", "K12D11" = "-4", "K12D13" = "-2", "K12D15" = "0", "K18D47" = "-10", "K18D49" = "-8", "K18D51" = "-6", "K18D53" = "-4", "K18D55" = "-2", "K18D57" = "0", "K19D23" = "-10", "K19D25" = "-8", "K19D27" = "-6", "K19D29" = "-4", "K19D31" = "-2", "K19D33" = "0", "K20D75" = "-10", "K20D77" = "-8", "K20D79" = "-6", "K20D81" = "-4", "K20D83" = "-2", "K20D85" = "0")


Day_labels=c("K12D04" = "-11", "K12D06" = "-9", "K12D09" = "-6", "K12D11" = "-4", "K12D13" = "-2", "K12D15" = "0", "K18D47" = "-10", "K18D49" = "-8", "K18D51" = "-6", "K18D53" = "-4", "K18D55" = "-2", "K18D57" = "0", "K19D23" = "-10", "K19D25" = "-8", "K19D27" = "-6", "K19D29" = "-4", "K19D31" = "-2", "K19D33" = "0", "K20D75" = "-10", "K20D77" = "-8", "K20D79" = "-6", "K20D81" = "-4", "K20D83" = "-2", "K20D85" = "0")





################################### SHOULDNT NEED TO MESS WITH THIS SECTION ########################
                                  


biobake_tax_input_file <- paste(biobake_input_folder,"/metaphlan/merged/metaphlan_taxonomic_profiles.tsv", sep = "")


Metaphlan_to_VALENCIA_output <- paste( Metaphlan_analysis_output_folder, Metaphlan_analysis_outputfilename, sep = "")

Metaphlan_CST_vis_input <- paste(Valencia_output_file_destination, "/", Valencia_output, sep = "")
                                      
                                      
################### Actual Code ##########################################
                                      

 ###use paste to save variables for file names


taxo_data <- read.csv( biobake_tax_input_file, sep="\t", stringsAsFactors=FALSE) %>% 
  rename(taxonomy = X..taxonomy ) 




##Fixing Species Names 

taxo_data <- filter(taxo_data, grepl('s__', taxonomy)) 

taxo_data$taxonomy <- gsub(".*s__","",taxo_data$taxonomy) 
taxo_data_rotate <- t(taxo_data)
taxo_data_rotate <- as.data.frame(taxo_data_rotate) %>%
  rownames_to_column(var = "sampleID") %>%
  row_to_names(row_number = 1) %>% 
  filter(!grepl("R2_", taxonomy)) %>% 
  filter(!grepl("Mock", taxonomy)) %>% 
  filter(!grepl("commun", taxonomy)) %>% 
  filter(!grepl("QH", taxonomy))






taxo_data_rotate$sampleID <- taxo_data_rotate$taxonomy

taxo_data_rotate <- taxo_data_rotate %>%
  dplyr::select(-taxonomy) %>%
  mutate( sampleID =  str_sub(sampleID, 1,6))

taxo_data_rotate$sampleID[taxo_data_rotate$sampleID == "K12D4_"] <-  "K12D04"
taxo_data_rotate$sampleID[taxo_data_rotate$sampleID == "K12D6_"] <-  "K12D06"
taxo_data_rotate$sampleID[taxo_data_rotate$sampleID == "K12D9_"] <-  "K12D09"

taxo_data_rotate <- taxo_data_rotate[c(3,4,5,1,2,6:nrow(taxo_data_rotate)), ]

#taxo_data_rotate$read_count <- rowSums(taxo_data_rotate[ ,2:90], na.rm = TRUE)
rownames(taxo_data_rotate) <-  taxo_data_rotate$sampleID




#### Tidy Data, calculating Relative Abundance, and Subsetting Taxa by Abundance Rank

taxo_data_rotate1 <- taxo_data_rotate %>%
  pivot_longer(cols = !sampleID, names_to = "Taxonomy", values_to = "Relative_Abundance") %>%
  mutate( Patient = str_sub( sampleID, 1,3 ))
#group_by(sampleID) %>% 
#mutate(sample_sum = sum(as.numeric(Relative_Abundance)))

taxo_data_rotate2 <- taxo_data_rotate1 %>% 
  
  group_by(Taxonomy) %>% 
  summarize( tax_abun = sum(as.numeric(Relative_Abundance), .drop=FALSE)) %>%
  arrange(desc(tax_abun)) %>% 
  mutate(rank = row_number()) %>% 
  inner_join(taxo_data_rotate1) %>% 
  
  
  mutate(Taxonomy = if_else(rank > 10, "Other", Taxonomy )) %>%
  group_by(sampleID) %>% 
  
  mutate(sample_sum = sum(as.numeric(Relative_Abundance))) %>%
  ungroup() 


###Color Stuff, dont think this is working yet 

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color_hue(length(levels(taxo_data_rotate2$Taxonomy)))

cols[levels(taxo_data_rotate2$Taxonomy)=="Other"]="grey"





################ Tax Abundance Bar Plot

ggplot(taxo_data_rotate2, aes(x = sampleID, y = as.numeric(Relative_Abundance), fill = Taxonomy, color=I("black"))) +
  geom_col() +
  theme_classic() + facet_grid(~Patient, scales = "free") +
  labs(x = "Days Prior to iBV", y = "Relative Abundance") +
  scale_x_discrete(labels=xlabels2)+
  theme(legend.text = element_text(face = "italic"))+
  scale_y_continuous(expand=c(0,0), limits=c(0,100.05))




################# Tax Heatmap

ggplot(taxo_data_rotate1, aes(x = sampleID, y = Taxonomy) )+
  geom_raster(aes(fill = as.numeric(Relative_Abundance))) +
  theme_classic() + facet_grid(~Patient, scales = "free") +
  theme(axis.text.y.left = element_text(face = "italic"))+
  scale_fill_gradient(name = NULL) +
  scale_x_discrete(labels=xlabels2) 



################ Metaphlan To Valencia ############################



taxo_data_rotate_met <- t(taxo_data)
taxo_data_rotate_met <- as.data.frame(taxo_data_rotate_met) %>%
  row_to_names(row_number = 1) 

taxo_data_rotate_met1 <- as.data.frame(sapply(taxo_data_rotate_met, as.numeric)) %>% 
  mutate( read_count = rowSums(.)) %>% 
  mutate_all(~100*(.))

rownames(taxo_data_rotate_met1) <- rownames(taxo_data_rotate_met)
taxo_data_rotate_met2 <- taxo_data_rotate_met1 %>% rownames_to_column(var = "sampleID") %>% 
  mutate( sampleID =  str_sub(sampleID, 1,6)) %>% 
  dplyr::select(sampleID, read_count, everything(), ) 


### Executing Terminal Commands To Run VALENCIA
session <- ssh_connect("jlamm1@tigerfish")
write.csv(taxo_data_rotate_met2, Metaphlan_to_VALENCIA_output, row.names = FALSE)
scp_upload(session, Metaphlan_to_VALENCIA_output )
ssh_exec_wait(session, command = c( paste( 'cp', Metaphlan_analysis_outputfilename,  'VALENCIA'),
                                    'cd VALENCIA',
                                    'conda activate VIRGO',
                                    paste('srun python3 ./Valencia.py -ref CST_centroids_012920.csv -i', Metaphlan_analysis_outputfilename, '-o', Valencia_output),
                                    paste('cp ', Valencia_output,  ' ../', sep="" )
                                    )
                                     )

        
                  

scp_download(session, Valencia_output )
system(paste('mv', Valencia_output, Valencia_output_file_destination ))
ssh_disconnect(session)

################################################# CST Analysis ##################################



Metaphlan_CST <- read.csv(Metaphlan_CST_vis_input) %>% 
  dplyr::select(sampleID, score, CST)


taxo_meta_dat <- taxo_data_rotate %>%
  inner_join(Metaphlan_CST, by = "sampleID") %>% 
  mutate(Patient = str_sub(sampleID, 1,3))


Day_labels <- Day_labels %>% as.data.frame() %>% rownames_to_column(var = "sampleID") %>% rename(Day = ".")


######### Creating Dist Matrix########

CST_mat <- taxo_data_rotate %>% dplyr::select(-sampleID)

CST_mat <- sapply(CST_mat, as.numeric)

CST_mat <- as.matrix(CST_mat)

set.seed(1)
dist <- vegdist(CST_mat, method = "bray")
nmds <- metaMDS(dist)

goodness(nmds)
stressplot(nmds)
Metaphlan_nmds <- scores(nmds, display="site") %>%
  as_tibble(rownames = "sampleID") 
Metaphlan_nmds$sampleID <- taxo_meta_dat$sampleID

metadata_nmds <- inner_join(Metaphlan_nmds, taxo_meta_dat, by= "sampleID")  %>% 
  inner_join(Day_labels, by = "sampleID")



metadata_nmds %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=Patient, label = Day ))+
  geom_label_repel()+
  geom_point( aes(shape =CST, size = score))+
  scale_shape_manual(values=c(6,16,18,17))+
  theme_classic()




