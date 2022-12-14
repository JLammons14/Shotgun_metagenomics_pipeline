library(tidyverse)
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
library(viridis)
library(circlize)
library(ComplexHeatmap)
library(hrbrthemes)
library(fs)
library(data.table)
library(Rpython)
library(jcolors)
library(ggnewscale)





###Instructions: VIRGO and VALENCIA Folders should be your cluster profile home directory, if not we can route to it in the future
                ### I run these programs while in a conda env on the cluster. 


###Automating the file Transfer
#################### THIS IS YOUR INPUT/OUTPUT SECTION #####################

Project_name <- "K_shotgun_"
output_folder <- "~/VIRGO_Valencia_Results_K_shotgun"

Raw_reads_location <- "../../media/scratch/john/shotgun_data"

Worflows_output_folder_name <- "biobakery_2"

username <- jlamm1

biobake_input_folder <- "~/biobakery_2"
Nugent_meta_data <- read.csv("~/Documents/Kshotgun_nuget_metadata.csv",sep = "\t", header = TRUE)

Cluster_conda_env <- "VIRGO"
#####################################Less important Inputs 

Metaphlan_analysis_outputfilename <- paste(Project_name,"metaphlan.csv", sep="")

Metaphlan_analysis_output_folder <- paste(output_folder, "/VIRGO_Valencia_Input/", sep="")


Valencia_output <- paste(Project_name,"metaphlan_CST.csv", sep="")

Valencia_output_file_destination <- paste(output_folder, "/Valencia", sep="")




VIRGO_local_folder <- paste(output_folder,"/VIRGO_K_output", sep="")



VIRGO_to_Valencia_folder <- paste(output_folder,"/VIRGO_Valencia_Input/", sep="")

Valencia_analysis_of_VIRGO_tax <-paste(Project_name,"VIRGO_CST", sep="")


 


xlabels3=c("K12D04" = "-11", "K12D06" = "-9", "K12D09" = "-6", "K12D11" = "-4", "K12D13" = "-2", "K12D15" = "0", "K18D47" = "-10", "K18D49" = "-8", "K18D51" = "-6", "K18D53" = "-4", "K18D55" = "-2", "K18D57" = "0", "K19D23" = "-10", "K19D25" = "-8", "K19D27" = "-6", "K19D29" = "-4", "K19D31" = "-2", "K19D33" = "0", "K20D75" = "-10", "K20D77" = "-8", "K20D79" = "-6", "K20D81" = "-4", "K20D83" = "-2", "K20D85" = "0")

xlabels2=c("K12D04" = "-11", "K12D06" = "-9", "K12D09" = "-6", "K12D11" = "-4", "K12D13" = "-2", "K12D15" = "0", "K18D47" = "-10", "K18D49" = "-8", "K18D51" = "-6", "K18D53" = "-4", "K18D55" = "-2", "K18D57" = "0", "K19D23" = "-10", "K19D25" = "-8", "K19D27" = "-6", "K19D29" = "-4", "K19D31" = "-2", "K19D33" = "0", "K20D75" = "-10", "K20D77" = "-8", "K20D79" = "-6", "K20D81" = "-4", "K20D83" = "-2", "K20D85" = "0")


Day_labels=c("K12D04" = "-11", "K12D06" = "-9", "K12D09" = "-6", "K12D11" = "-4", "K12D13" = "-2", "K12D15" = "0", "K18D47" = "-10", "K18D49" = "-8", "K18D51" = "-6", "K18D53" = "-4", "K18D55" = "-2", "K18D57" = "0", "K19D23" = "-10", "K19D25" = "-8", "K19D27" = "-6", "K19D29" = "-4", "K19D31" = "-2", "K19D33" = "0", "K20D75" = "-10", "K20D77" = "-8", "K20D79" = "-6", "K20D81" = "-4", "K20D83" = "-2", "K20D85" = "0")

Day_labels_verbose=c("K12D04" = "K12D-11", "K12D06" = "K12D-9", "K12D09" = "K12D-6", "K12D11" = "K12D-4", "K12D13" = "K12D-2", "K12D15" = "K12D0", "K18D47" = "K18D-10", "K18D49" = "K18D-8", "K18D51" = "K18D-6", "K18D53" = "K18D-4", "K18D55" = "K18D-2", "K18D57" = "K18D0", "K19D23" = "K19D-10", "K19D25" = "K19D-8", "K19D27" = "K19D-6", "K19D29" = "K19D-4", "K19D31" = "K19D-2", "K19D33" = "K19D0", "K20D75" = "K20D-10", "K20D77" = "K20D-8", "K20D79" = "K20D-6", "K20D81" = "K20D-4", "K20D83" = "K20D-2", "K20D85" = "K20D0")





################################### SHOULDNT NEED TO MESS WITH THIS SECTION ########################
                                  


biobake_tax_input_file <- paste(biobake_input_folder,"/metaphlan/merged/metaphlan_taxonomic_profiles.tsv", sep = "")


Metaphlan_to_VALENCIA_output <- paste( Metaphlan_analysis_output_folder, Metaphlan_analysis_outputfilename, sep = "")

Metaphlan_CST_vis_input <- paste(Valencia_output_file_destination, "/", Valencia_output, sep = "")

Human_heatmap_vis_input <- paste(biobake_input_folder, "/humann/merged/pathabundance_relab.tsv", sep ="")

Kneaddata_outputfiles <- paste(biobake_input_folder, "/kneaddata/main", sep="") 

VIRGO_vis_input_total_abun <- paste(VIRGO_local_folder, "/temp_mapping/summary.Abundance.txt", sep = "")

VIRGO_vis_input_Rel_abund <- paste(VIRGO_local_folder, "/temp_mapping/summary.Percentage.txt", sep = "")


VIRGO_to_valencia_file <- paste(Project_name, "VIRGO_to_Valencia.csv", sep = "")

VIRGO_local_folder_for_path_analysis <- paste(VIRGO_local_folder, "/temp_mapping", sep="")


wd <- getwd()
local_user <- system(" id -un", intern = TRUE)
machine_name <- system("hostname", intern = TRUE)
cluster_name <- "tigerfish"
                                      
################### Actual Code ##########################################

 ###Biobakery Workflows Command ### You may want to run this step just in the terminal, but heres the generic code to run workflows 
#system(paste("biobakery_workflows wmgx --input", Raw_reads_location, "--output", Worflows_output_folder_name,  "--local-job 5 --threads 8  --bypass-strain-profiling"))

                                      


################### Tax Abundance Bar plot #############

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


taxo_data_rotate <- taxo_data_rotate[c(4,5,6,1,2,3,7:nrow(taxo_data_rotate)), ]

#taxo_data_rotate$read_count <- rowSums(taxo_data_rotate[ ,2:90], na.rm = TRUE)
rownames(taxo_data_rotate) <-  taxo_data_rotate$sampleID




#### Tidy Data, calculating Relative Abundance, and Subsetting Taxa by Abundance Rank

taxo_data_rotate1 <- taxo_data_rotate %>% 
  pivot_longer(cols = !sampleID, names_to = "Taxonomy", values_to = "Relative_Abundance") %>%
  mutate( Patient = str_sub( sampleID, 1,3 )) %>% 
  mutate(Taxonomy = replace(Taxonomy, Taxonomy == "Atopobium_vaginae", "Fannyhessea_vaginae")) %>% 
  mutate(category = "2")

taxo_data_rotate$sampleID[taxo_data_rotate$sampleID == "K12D4_"] <-  "K12D04"
#group_by(sampleID) %>% 
#mutate(sample_sum = sum(as.numeric(Relative_Abundance)))


taxo_data_rotate2 <- taxo_data_rotate1 %>% 
  
  group_by(Taxonomy) %>% 
  summarize( tax_abun = sum(as.numeric(Relative_Abundance), .drop=FALSE)) %>%
  arrange(desc(tax_abun)) %>% 
  mutate(rank = row_number()) %>% 
  inner_join(taxo_data_rotate1) %>% 
  
  
  mutate(Taxonomy = if_else(rank > 12, "Other", Taxonomy )) %>%
  group_by(sampleID) %>% 

  
  mutate(sample_sum = sum(as.numeric(Relative_Abundance))) %>%
  ungroup() %>% 
  mutate(Taxonomy = fct_reorder(Taxonomy, rank, .desc = TRUE)) %>% 
  filter(sampleID !="K12D11")


taxo_data_rotate2_adjusted <- taxo_data_rotate2 %>% 
  mutate(Relative_Abundance = as.numeric(Relative_Abundance)/100)


###Color Stuff, dont think this is working yet 

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color_hue(length(levels(taxo_data_rotate2$Taxonomy)))

cols[levels(taxo_data_rotate2$Taxonomy)=="Other"]="grey"





################ Tax Abundance Bar Plot
stack_barplot <- function(tax_df, sampleID1, values1, tax_names, facet_var, sample_labels, ylab) {

p <- ggplot(tax_df, aes(x = sampleID1, y = as.numeric(values1), fill = tax_names, color=I("black"))) +
  geom_col() +
  theme_classic() + facet_grid(facet_var, scales = "free") +
  labs(x = "Days Prior to iBV", y = ylab) +
  scale_x_discrete(labels=sample_labels)+
  theme(legend.text = element_text(face = "italic"))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_discrete(name = "Taxonomy")
  return(p)
}

stack_barplot(taxo_data_rotate2, taxo_data_rotate2$sampleID, taxo_data_rotate2$Relative_Abundance, taxo_data_rotate2$Taxonomy, ~Patient, xlabels2, "Relative Abundance")




################# Tax Heatmap
Tax_heatmap <- function(tax_df, sampleID, Tax, values, groups1, sample_labels) {

ggplot(tax_df, aes(x = sampleID, y = Tax) )+
  geom_raster(aes(fill = log10(as.numeric(values)+0))) +
  theme_classic() + facet_grid(groups1, scales = "free") +
  theme(axis.text.y.left = element_text(face = "italic"))+
   # scale_fill_gradient(breaks = c(-7.5, -6, -4, -2, -1), labels = c("0",".0001", ".01", "1", "10"), na.value="black", low = "darkblue", high = "lightblue", name = "Bacteria abundance (%)", guide = guide_colourbar(direction = "horizontal", order = 2, title.position = "top")) +
  scale_x_discrete(labels=sample_labels) +
    labs(x = "Days Prior to iBV", y = "Taxonomy")
  
}

Tax_heatmap(taxo_data_rotate1, sampleID, factor(taxo_data_rotate1$Taxonomy), taxo_data_rotate1$Relative_Abundance, ~Patient, xlabels2)
 Tax_heatmap(taxo_data_rotate2_adjusted, sampleID, factor(taxo_data_rotate2_adjusted$Taxonomy), taxo_data_rotate2_adjusted$Relative_Abundance, ~Patient, xlabels2)



################## Custom Heatmap for addining in Nuget Scores ################################


Nugent_meta_data <- Nugent_meta_data %>% mutate(Patient = str_sub(sampleID, 1,3)) %>% 
  filter(sampleID != "K12D11")
 
 biobake_tax_heatmap <- taxo_data_rotate2_adjusted  %>% group_by(Taxonomy, sampleID) %>% 
   mutate(rel_abund_heatmap = sum(as.numeric(Relative_Abundance))) %>% 
   ungroup() %>% 
   mutate(category = 2) %>% 
   filter(sampleID != "K12D11") %>% 
   mutate(Taxonomy = fct_reorder(Taxonomy, rank, .desc = TRUE))
 
  
 ggplot() +
  geom_raster(data = biobake_tax_heatmap, aes(x = sampleID, y = Taxonomy, fill = as.numeric(log10(rel_abund_heatmap))))  +
   theme(
    strip.text.y = element_blank()) +
facet_grid( rows = vars(category),
            cols = vars(Patient), scales = "free", space = "free_y") +
  
  theme(axis.text.y.left = element_text(face = "italic"))+
  scale_x_discrete(labels=xlabels2) +
  labs(x = "Days Prior to iBV", y = "Taxonomy") +
   
   scale_fill_gradient(breaks = c( -5, -3,  -2, -1), 
                       labels = c(".0001",  ".01", "1", "10" ),
                       na.value="black",
                       low = "darkblue", high = "lightblue", 
                       name = "Bacteria abundance (%)", 
                       guide = guide_colourbar(direction = "horizontal", 
                                               #order = 2,
                                               title.position = "top")
                       ) +
   
   
   new_scale_fill() +
   
   geom_raster(data = Nugent_meta_data, aes(x = sampleID, y = Score, fill = Nuget_Category ))  +
  
   scale_fill_manual(values = c(Low = "green", Mid = "yellow",  High = "red"), name = "Nuget Score") +
   geom_text(data = Nugent_meta_data, aes(x = sampleID, y = Score, label = Nuget_Score), size = 4)+
   geom_text(data = Nugent_meta_data, aes(x = sampleID, y = Score, label = Nuget_Score), size = 4)+
  #  geom_text(data = Nugent_meta_data, aes(x = sampleID, y = Score, label = Nuget_Score), size = 4)+
   guides(fill = "none")


 # # guides(fill =
 #           guide_legend(
 #             title.theme = element_text(
 #               size = 10), label.hjust = .05))

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
session <- ssh_connect(paste(local_user, "@", cluster_name, sep = ""))
write.csv(taxo_data_rotate_met2, Metaphlan_to_VALENCIA_output, row.names = FALSE)
scp_upload(session, Metaphlan_to_VALENCIA_output, to ="VALENCIA" )
ssh_exec_wait(session, command = c( 
                                    'cd VALENCIA',
                                    paste('conda activate',  Cluster_conda_env),
                                    paste('srun python3 ./Valencia.py -ref CST_centroids_012920.csv -i', Metaphlan_analysis_outputfilename, '-o', Valencia_output)
                                  
                                    )
                                     )

        
                  

scp_download(session, paste('VALENCIA/', Valencia_output, sep="" ), to = paste( Valencia_output_file_destination, sep = ""))
ssh_disconnect(session)

################################################# CST Analysis ##################################



Metaphlan_CST <- read.csv(Metaphlan_CST_vis_input) %>% 
  dplyr::select(sampleID, score, CST)


taxo_meta_dat <- taxo_data_rotate %>%
  inner_join(Metaphlan_CST, by = "sampleID") %>% 
  mutate(Patient = str_sub(sampleID, 1,3))


Day_labels <- Day_labels %>% as.data.frame() %>% rownames_to_column(var = "sampleID") %>% rename(Day = ".")


######### Creating Dist Matrix######
  


CST_mat <- taxo_data_rotate %>% dplyr::select(-sampleID) 

CST_mat2 <- CST_mat[-c(1),]

CST_mat <- sapply(CST_mat2, as.numeric)

CST_mat1 <- as.matrix(CST_mat)

Metaphlan_CST_add_count <- as.data.frame(CST_mat) %>% 
  mutate(count = rowSums(.)) %>% 
  rownames_to_column(var="sampleID") %>% 
  dplyr::select(sampleID, count)

set.seed(1)
dist <- vegdist(CST_mat1, method = "bray")
nmds <- metaMDS(dist)


goodness(nmds)
stressplot(nmds)
Metaphlan_nmds <- scores(nmds, display="site") %>%
  as_tibble(rownames = "sampleID") 
Metaphlan_nmds$sampleID <- rownames(CST_mat2)

metadata_nmds <- left_join(Metaphlan_nmds, taxo_meta_dat, by= "sampleID")  %>% 
  left_join(Day_labels, by = "sampleID")

metadata_nmds %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=Patient, label = Day ))+
  geom_label_repel(key_glyph = "rect")+
  geom_point( aes(shape =CST, size = 10))+
  scale_shape_manual(values=c("\u29EB","\u25FC","\u25B2","\u25CF"))+
  coord_flip() +
  theme_classic() +
   scale_color_manual(values = c("#a9d1e3", "#8c33e4", "#00008B", "#000000"))

################### Ordination Alternative 
tax_dist <- dist %>% as.matrix() %>% as.data.frame()
rownames(tax_dist) <- rownames(CST_mat2)
tax_dist1 <- tax_dist %>% mutate(sampleID = rownames(tax_dist)) %>% 
                          
            pivot_longer(-sampleID, values_to = "distance", names_to = "names" ) %>% 
            left_join(Metaphlan_CST, by = "sampleID") %>% 
            mutate(Patient = str_sub(sampleID, 1,3))
  
  
  
  

############### Vis Humann Pathway Analysis ##################


path <- read.csv(Human_heatmap_vis_input, sep="\t") %>%
  
  #path <- read.csv("../Downloads/pathabundance.tsv", sep="\t")   %>%
  rename(Pathway = X..Pathway) 

path_short <- path %>% 
  filter( !grepl('g__', Pathway)) %>% 
  dplyr::select(-contains("QH"), 
                -contains("community"),
                -contains("R2_")) 
  

path_shot_test1 <- path_short %>%
  mutate(data.frame(str_split_fixed(path_short$Pathway, ":", 2))) %>% 
  rename(Pathway_Name = X2) %>% 
  rename(Pathway_symbol = X1) %>% 
  group_by(Pathway_Name) %>% 
  mutate(path_relabund_sum = rowSums(across(where(is.numeric)))) 

path_shorter <- path_shot_test1 %>%
  #filter( !grepl('unclassified', Pathway_Name)) %>% 
  column_to_rownames( var = "Pathway_Name") %>%
  dplyr::select(-Pathway) %>%
  arrange(desc(path_relabund_sum)) %>%
  mutate(rank = row_number() ) %>%
  filter(rank <= 50) %>%
  dplyr::select (-path_relabund_sum, -rank, -Pathway_symbol)# %>% 
 #dplyr::select(K12D4_S25_L003_R1_001_Abundance, K12D6_S26_L003_R1_001_Abundance, K12D9_S27_L003_R1_001_Abundance, everything() )

path_rotate <-  t(path_shorter) %>%
  as.data.frame() %>%
  rownames_to_column("sampleID") %>%  
  mutate(labs = xlabels2)
#%>%
#pivot_longer( -sampleID, names_to = "Pathway", values_to = "Relative_Abundance")  


path_rotate$sampleID <- gsub("_.*","", path_rotate$sampleID)


rownames(path_rotate) <- path_rotate$sampleID

path_rotate1 <- path_rotate %>% dplyr::select( -sampleID, -labs) 

CST_meta <- metadata_nmds %>% as.data.frame() %>%
  dplyr::select(sampleID, CST) %>%
  arrange(sampleID) %>%
  mutate(ID = Day_labels_verbose) %>% 
  dplyr::select(-sampleID) %>% 
  rename(CST = CST)

rownames(CST_meta) <- CST_meta$ID
CST_meta<- CST_meta %>% dplyr::select(-ID) %>% 
  rownames_to_column(var="sampleID")

path_rotate1 <-t(path_rotate1)

colnames(path_rotate1) <- Day_labels_verbose

#CST_meta <- t(CST_meta)

paletteLength <- 50

mycolor <- viridis::viridis
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

path_rotate1 <- as.matrix(path_rotate1)

zscores <- scale_rows(path_rotate1)
ha_c = HeatmapAnnotation(CST = CST_meta$CST,  
                         col = list(CST = c("I"="#089392", "II" = "#D95F02", "III" = "#82C782", "IV-B" = "#CF597E", "V" = "#EACF87"))) 
   
                        
ht1 = Heatmap(zscores, name = " ",
              col = viridis(50),
              top_annotation = ha_c,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              row_names_max_width = unit(10, "cm"),
              heatmap_legend_param = list(title = " ", legend_height = unit(4, "cm")),
              column_split = CST_meta$CST
              ) 
draw(ht1 )
pheatmap( mat = zscores, 
         
         annotation_col = CST_meta, 
         color = viridis(50), 
       
         #legend=F, 
          #legend_labels = c("Z-Score", "CST"),
         #cluster_cols = FALSE,
         #gaps_col = c(6,12,18),
         border_color=NA
         
         
         
) 
draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))

################################# Completed Biobakery Portion of Analysis ##################################
################################# starting the VIRGO Section of Analysis ##################################
#################################                                       ##################################

### To Run the VIRGO portion of the analysis, I formerly utilized multiline editing 
### in Atom to create the input file. I want to automate this step as well, but it will take a bit 
## so for now I'll skip that step. 


session <- ssh_connect(paste(local_user, "@", cluster_name, sep = ""))


ssh_exec_wait(session, command = c( 
              
              'mkdir VIRGO/4_run_VIRGO/kneaddata_output'
              ))

### Moving cleaned files into a seperate directory so they can be transfered

system(paste("mkdir ", biobake_input_folder, '/kneaddata/main/clean_forward_reads', sep=""))
system(paste("mv ", biobake_input_folder, "/kneaddata/main/*R1_001.fastq ", biobake_input_folder, '/kneaddata/main/clean_forward_reads' , sep=""))
scp_upload(session, paste(Kneaddata_outputfiles,"/clean_forward_reads", sep = ""), to = "VIRGO/4_run_VIRGO/kneaddata_output")


###Cleaning up garbaged directory I made

system(paste("mv ", biobake_input_folder, "/kneaddata/main/clean_forward_reads/*R1_001.fastq ", biobake_input_folder, '/kneaddata/main' , sep=""))

system(paste("rm -r ", biobake_input_folder, '/kneaddata/main/clean_forward_reads' , sep="" ))
#ssh_exec_wait(session, 
ssh_exec_wait(session, command = c("rm -r VIRGO/4_run_VIRGO/temp_mapping"))
##Moving into VIRGO to RUN the program
#"sbatch ./runMapping.step1.sh -r ../shotgun_cleaned_files/K12D13_S29_L003_R1_001_kneaddata_paired_1.fastq -p K12D13  -d /mnt/beegfs/home/jlamm1/VIRGO")
ssh_exec_wait(session, command = c( "cd VIRGO/4_run_VIRGO/",
                          "ls kneaddata_output/clean_forward_reads > kneaddata_output_files.txt",
                           paste('conda activate', Cluster_conda_env),
                          "for file in $(cat kneaddata_output_files.txt); do sbatch ./runMapping.step1.sh -r kneaddata_output/clean_forward_reads/$file -p ${file::6}  -d ~/VIRGO ;  done"
))          

 

###paste("scp -r temp_mapping ", local_user, "@", machine_name, ":",VIRGO_output, sep = "")


### This is to prevent movement of files before they finish running

file_num <- 0
file_num_tot <- 1

while(file_num[1] != file_num_tot[1] ) {
   file_num <- capture.output(ssh_exec_wait(session, command = c(" num=$(wc -l < VIRGO/4_run_VIRGO/kneaddata_output_files.txt)", 
                                        "num_mult=$((($num)))",
                                          'echo  "$num_mult"') )) 
   file_num_tot <- capture.output( ssh_exec_wait(session, command = "ls -lq VIRGO/4_run_VIRGO/temp_mapping  | wc -l"  )) 
   Sys.sleep(300)}   



ssh_exec_wait(session, command = c( "cd VIRGO/4_run_VIRGO/",
                                    paste('conda activate', Cluster_conda_env),
                                    "sbatch ./runMapping.step2.sh -p temp_mapping -d ~/VIRGO"))


file_num <- 0
file_num_tot <- 1

while( file_num[1] != file_num_tot[1]) {
  file_num <- capture.output(ssh_exec_wait(session, command = c(" num=$(wc -l < VIRGO/4_run_VIRGO/kneaddata_output_files.txt)", 
                                     "num_mult=$((($num * 11)+6))",
                                     'echo  "$num_mult"') ) )
  
  file_num_tot <- capture.output(ssh_exec_wait(session, command = "ls -lq VIRGO/4_run_VIRGO/temp_mapping  | wc -l"  ))
  Sys.sleep(60)
  }   




  scp_download(session, "VIRGO/4_run_VIRGO/temp_mapping",  to = VIRGO_local_folder)
  ssh_exec_wait(session, command = "rm -r VIRGO/4_run_VIRGO/temp_mapping")
ssh_disconnect(session)



                
############################# WE did it! VIRGO complete!!!!!!!!!! Now plot abundance ######################

 ###### VIRGO Analysis Absolute Abundance barplot 

### Parts fo this were the first  code I wrote for the project, I came back and 
## cleaned some of it to make certain improvements,so the differences are suddle, but 
## see if you can notice the improvements I've made. 

VIRGO_taxo <- read.delim(VIRGO_vis_input_total_abun) %>% 
  filter(!grepl("Mock", PC)) %>% 

  filter(!grepl("Tth", PC))

VIRGO_taxo$PC[VIRGO_taxo$PC == "K12D4"] <-  "K12D04"
VIRGO_taxo$PC[VIRGO_taxo$PC == "K12D6"] <-  "K12D06"
VIRGO_taxo$PC[VIRGO_taxo$PC == "K12D9"] <-  "K12D09"
VIRGO_taxo <- VIRGO_taxo[c(3,4,5,1,2,6:nrow(VIRGO_taxo)), ]

#VIRGO_taxo$PC <- gsub("\\_.*", "", VIRGO_taxo$PC)

VIRGO_taxo1 <- VIRGO_taxo 

VIRGO_taxo1$Sample <- VIRGO_taxo1$PC
VIRGO_taxo1$Patient <- str_sub(VIRGO_taxo1$Sample, 1, 3)

names <- str(VIRGO_taxo1[1,])
VIRGO_taxo_longer <- VIRGO_taxo1 %>%
  pivot_longer(-c(PC,Sample,Patient), names_to = "Taxonomy", values_to = "Abundance") 

VIRGO_Combo <- VIRGO_taxo_longer %>% 
  group_by(Taxonomy) %>%
  summarize( tax_abun = sum(as.numeric(Abundance), .drop=FALSE)) %>%
  arrange(desc(tax_abun)) %>% 
  mutate(rank = row_number()) %>% 
  inner_join(VIRGO_taxo_longer) %>%
  mutate(Taxonomy = if_else(rank > 12, "Other", Taxonomy ))# %>%

VIRGO_Combo_vis <- dplyr::select(VIRGO_Combo, Sample, Abundance, Taxonomy, Patient)

VIRGO_xlabels2=c("K12D04" = "-11", "K12D06" = "-9", "K12D09" = "-6", "K12D13" = "-2", "K12D15" = "0", "K18D47" = "-10", "K18D49" = "-8", "K18D51" = "-6", "K18D53" = "-4", "K18D55" = "-2", "K18D57" = "0", "K19D23" = "-10", "K19D25" = "-8", "K19D27" = "-6", "K19D29" = "-4", "K19D31" = "-2", "K19D33" = "0", "K20D75" = "-10", "K20D77" = "-8", "K20D79" = "-6", "K20D81" = "-4", "K20D83" = "-2", "K20D85" = "0")

stack_barplot(VIRGO_Combo_vis, VIRGO_Combo_vis$Sample, VIRGO_Combo_vis$Abundance, VIRGO_Combo_vis$Taxonomy, ~Patient, VIRGO_xlabels2, "Abundance")


###################################################### VIRGO Bar Plot Rel Abund ##########


VIRGO_taxo_rel <- read.delim(VIRGO_vis_input_Rel_abund) %>% 
  filter(!grepl("Mock", PC)) %>% 
filter(!grepl("Tth", PC))


VIRGO_taxo_rel$PC[VIRGO_taxo_rel$PC == "K12D4"] <-  "K12D04"
VIRGO_taxo_rel$PC[VIRGO_taxo_rel$PC == "K12D6"] <-  "K12D06"
VIRGO_taxo_rel$PC[VIRGO_taxo_rel$PC== "K12D9"] <-  "K12D09"
VIRGO_taxo_rel <- VIRGO_taxo_rel[c(4,5,6,1,2,3,7:nrow(VIRGO_taxo_rel)), ]


VIRGO_taxo_rel$Sample <- VIRGO_taxo_rel$PC
VIRGO_taxo_rel$Patient <- str_sub(VIRGO_taxo_rel$Sample, 1, 3)

#VIRGO_taxo_rel$Sample <- gsub("\\_.*", "", VIRGO_taxo_rel$Sample)
names <- str(VIRGO_taxo_rel[1,])
VIRGO_taxo_longer_rel <- VIRGO_taxo_rel %>%
  pivot_longer(-c(PC,Sample,Patient), names_to = "Taxonomy", values_to = "Relative_Abundance") %>% 
  mutate(Taxonomy = replace(Taxonomy, Taxonomy == "Atopobium_vaginae", "Fannyhessea_vaginae" ))

VIRGO_Combo_rel <- VIRGO_taxo_longer_rel %>% 
  group_by(Taxonomy) %>%
  summarize( tax_abun_rel = sum(as.numeric(Relative_Abundance), .drop=FALSE)) %>%
  arrange(desc(tax_abun_rel)) %>% 
  mutate(rank = row_number()) %>% 
  inner_join(VIRGO_taxo_longer_rel) %>%
  mutate(Taxonomy = if_else(rank > 12, "Other", Taxonomy )) %>%
  mutate(Taxonomy = if_else(Taxonomy == "Coprobacillus_sp.", "Other", Taxonomy)) %>% 
  mutate(Taxonomy = fct_reorder(Taxonomy, rank, .desc = TRUE)) %>% 
    mutate( category = "2")

VIRGO_Combo_rel_heat <- VIRGO_Combo_rel %>% group_by(Taxonomy, Sample) %>% mutate(rel_abund_heatmap = sum(Relative_Abundance))

stack_barplot(VIRGO_Combo_rel, VIRGO_Combo_rel$Sample, VIRGO_Combo_rel$Relative_Abundance, VIRGO_Combo_rel$Taxonomy, ~Patient, VIRGO_xlabels2, "Relative Abundance")


Tax_heatmap(VIRGO_Combo_rel,VIRGO_Combo_rel$Sample ,VIRGO_Combo_rel$Taxonomy, VIRGO_Combo_rel$Relative_Abundance,~Patient, VIRGO_xlabels2  )



################################ CUSTOM NUGET HEATMAP ###########################

Nugent_meta_data <- Nugent_meta_data %>% mutate(Patient = str_sub(sampleID, 1,3)) %>% 
  filter(sampleID != "K12D11")

biobake_tax_heatmap <- taxo_data_rotate2_adjusted %>% group_by(Taxonomy, sampleID) %>% 
  mutate(rel_abund_heatmap = sum(as.numeric(Relative_Abundance))) %>% 
  mutate(Taxonomy = fct_reorder(Taxonomy, rank, .desc = TRUE)) %>%
  mutate(category = 2) %>% 
  filter(sampleID != "K12D11")

ggplot() +
  geom_raster(data = biobake_tax_heatmap, aes(x = sampleID, y = factor(Taxonomy), fill = log10(as.numeric(rel_abund_heatmap)+0)))  +
  scale_fill_gradient(low = "darkblue", high = "lightblue", name = "Relative Abundance")+
  theme(
    strip.text.y = element_blank()) +
  
  facet_grid( rows = vars(category),
              cols = vars(Patient), scales = "free", space = "free_y") +
  
  theme(axis.text.y.left = element_text(face = "italic"))+
  scale_x_discrete(labels=xlabels2) +
  labs(x = "Days Prior to iBV", y = "Taxonomy") +
  
  new_scale_fill() +
  
  geom_raster(data = Nugent_meta_data, aes(x = sampleID, y = Score, fill = Nuget_Category ))  +
  
  scale_fill_manual(values = c(Low = "green", Mid = "yellow",  High = "red"), name = "Nuget Score") +
  geom_text(data = Nugent_meta_data, aes(x = sampleID, y = Score, label = Nuget_Score), size = 4)+
  geom_text(data = Nugent_meta_data, aes(x = sampleID, y = Score, label = Nuget_Score), size = 4)+
  #  geom_text(data = Nugent_meta_data, aes(x = sampleID, y = Score, label = Nuget_Score), size = 4)+
  guides(fill = "none")





################################## TIME FOR VIRGO TO valencia ################################


###Table Prep
VIRGO_table_for_valencia <- VIRGO_taxo %>% 
mutate(read_count = rowSums(across(where(is.numeric)))) %>% 
rename(sampleID = PC) %>% 
dplyr::select(sampleID, read_count, everything(), )
write.csv(VIRGO_table_for_valencia, row.names = FALSE)


session <- ssh_connect(paste(local_user, "@", cluster_name, sep = ""))
write.csv(VIRGO_table_for_valencia, VIRGO_to_valencia_file, row.names = FALSE)
scp_upload(session, VIRGO_to_valencia_file, to ="VALENCIA" )
ssh_exec_wait(session, command = c( 
  'cd VALENCIA',
  paste('conda activate',  Cluster_conda_env),
  paste('srun python3 ./Valencia.py -ref CST_centroids_012920.csv -i', VIRGO_to_valencia_file, '-o', Valencia_analysis_of_VIRGO_tax)))




scp_download(session, paste('VALENCIA/', Valencia_analysis_of_VIRGO_tax, ".csv" , sep="" ), to = paste( Valencia_output_file_destination, sep = ""))
ssh_disconnect(session)


#############################ORD PLOT!! ###################################################



VIRGO_CST <- read.csv(paste(Valencia_output_file_destination,"/", Valencia_analysis_of_VIRGO_tax, ".csv", sep ="")) %>% 
  dplyr::select(sampleID, score, CST)



##creating metadata frame
VIRGO_taxo <- VIRGO_taxo %>% rename( sampleID = PC)


VIRGO_CST_meta <- VIRGO_taxo %>% 
  inner_join(VIRGO_CST, by = "sampleID") %>% 
  mutate(Patient = str_sub(sampleID, 1,3)) 



rownames(VIRGO_taxo) <- VIRGO_taxo$sampleID




VIRGO_xlabels2 <- xlabels2 %>% as.data.frame() %>% rownames_to_column(var = "sampleID") %>% rename(Day = ".")
VIRGO_xlabels3 <- xlabels3 %>% as.data.frame() %>% rownames_to_column(var = "sampleID") %>% rename(Day = ".")


##############make dist matrix #######################


VIRGO_CST_mat <- VIRGO_taxo %>% dplyr::select(-sampleID)

VIRGO_CST_mat1 <- sapply(VIRGO_CST_mat, as.numeric)

VIRGO_CST_mat1 <- as.matrix(VIRGO_CST_mat1)

VIRGO_CST_with_Count <- as.data.frame(VIRGO_CST_mat) %>% 
  mutate(count = rowSums(.)) %>% 
  rownames_to_column(var = "sampleID")


set.seed(1)
VIRGO_dist <- avgdist(VIRGO_CST_mat1, dmethod = "bray", sample=2000, iterations = 200)
VIRGO_nmds <- metaMDS(VIRGO_dist)

goodness(VIRGO_nmds)
stressplot(VIRGO_nmds)

VIRGO_nmds_df <- scores(VIRGO_nmds, display="site") %>%
  as_tibble(rownames = "sampleID") 

VIRGO_nmds_df$sampleID <- VIRGO_CST_meta$sampleID


VIRGO_metadata_nmds <- inner_join(VIRGO_nmds_df, VIRGO_CST_meta, by= "sampleID") %>%
  left_join(VIRGO_xlabels2, by = "sampleID") %>% 
  left_join(VIRGO_CST_with_Count)



VIRGO_metadata_nmds %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=Patient, label = Day ))+
  geom_label_repel(key_glyph = "rect")+
  geom_point( aes(shape =CST, size = score))+
  scale_shape_manual(values=c(6,16,18,17))+
 scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  guides(shape = guide_legend(order = 1),
         color  = guide_legend(order = 2),
          size = guide_legend(order = 0)) 
  #scale_color_jcolors(palette = "pal9")
  
  

####################################################LETS GO HEATMAP###############################
############################################################################################
#########################################################################################

setwd(paste(VIRGO_local_folder_for_path_analysis))

filenames_genelst <- list.files( pattern="*gene.lst.txt") 


VIRGO_gene_list = lapply(filenames_genelst, function(i){
  x=read_delim(i)
  x$sampleID = str_extract(i, "[^.]+")
  x
  
})

##### Loading in KEGG Gene ID's Basically

filenames_kegg.module.annotation <- list.files( pattern="*kegg.module.annotation.txt") 


VIRGO_kegg.module.annotation = lapply(filenames_kegg.module.annotation, function(i){
  x=read_delim(i)
  x$file = str_extract(i, "[^.]+")
  x
  
})

######## Loading in KEGG Pathway info 

filenames_kegg.pathway.annotation.txt <- list.files( pattern="*kegg.pathway.annotation.txt") 


VIRGO_kegg.pathway.annotation.txt = lapply(filenames_kegg.pathway.annotation.txt, function(i){
  x=read_delim(i)
  x$file = i
  x
  
})


##### Loading in Files for Protein Family info 


filenames_protein.family.annotation <- list.files( pattern="*proteinFamily.annotation.txt") 


VIRGO_protein.family.annotation = lapply(filenames_protein.family.annotation, function(i){
  x=read_delim(i)
  x$file = str_extract(i, "[^.]+")
  x
  
})

VIRGO_kegg.pathway.annotation.txt = do.call("rbind.data.frame", VIRGO_kegg.pathway.annotation.txt)

VIRGO_kegg.module.annotation = do.call("rbind.data.frame", VIRGO_kegg.module.annotation)

VIRGO_protein.family.annotation = do.call("rbind.data.frame", VIRGO_protein.family.annotation)




VIRGO_gene_info = do.call("rbind.data.frame", VIRGO_gene_list) %>%
  rename(Gene = gene)%>%
  filter(count > 200)  %>%
  left_join(  VIRGO_kegg.module.annotation, by = "Gene") %>%
  left_join( VIRGO_kegg.pathway.annotation.txt, by = "Gene") %>%
  left_join(VIRGO_protein.family.annotation, by = "Gene") %>% 
  filter(KEGG_module_number != "NA") 


VIRGO_kegg.pathway.annotation.txt = do.call("rbind.data.frame", VIRGO_kegg.pathway.annotation.txt)

VIRGO_kegg.module.annotation = do.call("rbind.data.frame", VIRGO_kegg.module.annotation)

VIRGO_protein.family.annotation = do.call("rbind.data.frame", VIRGO_protein.family.annotation)




VIRGO_gene_info = do.call("rbind.data.frame", VIRGO_gene_list) %>%
  rename(Gene = gene)%>%
  filter(count > 200)  %>%
  left_join(  VIRGO_kegg.module.annotation, by = "Gene") %>%
  left_join( VIRGO_kegg.pathway.annotation.txt, by = "Gene") %>%
  left_join(VIRGO_protein.family.annotation, by = "Gene") %>% 
  filter(KEGG_module_number != "NA") 
#filter(KEGG_module_number != "NA" |KEGG_pathwayID != "NA" ) %>%





VIRGO_gene_info$Patient <-  str_sub(VIRGO_gene_info$sampleID, 1, 3) 



VIRGO_gene_info <- VIRGO_gene_info %>% 
  dplyr::select(-file.x, -file.y, -file) %>%
  distinct() %>% 
  group_by(sampleID, KEGG_module_number) %>%
  filter(count > 100) %>% 
  mutate(pathway_count = sum(count)) %>%
  
  ungroup() %>% 
  group_by(sampleID) %>% 
  mutate(relabund = count/(sum(count))) %>%
  ungroup() 

VIRGO_gene_info$GO <-  str_sub(VIRGO_gene_info$GO, 3, -1)

VIRGO_sub <-  VIRGO_gene_info %>% 
  
  group_by(KEGG_module_number) %>% 
  mutate( total_path_count = sum(count)) %>% 
  ungroup() %>%  
  group_by(sampleID) %>% 
  summarise(rel_total_path_count = (total_path_count/(sum(count))), KEGG_module_number = KEGG_module_number, Gene = Gene) %>% 
  arrange(desc(rel_total_path_count)) %>% 
  mutate(rank = row_number()) %>%
  ungroup() %>% 
  inner_join(VIRGO_gene_info, by = c("sampleID" = "sampleID", "Gene" = "Gene")) %>%
  dplyr::select(Gene, count, rank, sampleID, KEGG_module_number.x,Annotation.x, GO, rel_total_path_count, count, KEGG_pathwayID, Annotation.y, relabund, Patient) %>% 
  group_by(sampleID, KEGG_module_number.x) %>% 
  mutate(sample_path_count= sum(count)) %>% 
  ungroup() %>% 
  group_by(sampleID) %>% 
  mutate(subject_sum_count = sum(count)) %>% 
  mutate(rel_sample_path_count = (sample_path_count/sum(count))) %>% 
  filter(rank < 30)

VIRGO_sub$GO <- str_split_fixed(VIRGO_sub$GO, ";", 2)[, 1]
VIRGO_sub$Subject <-  str_sub(VIRGO_sub$sampleID, 1, 3) 


VIRGO_sub$sampleID[VIRGO_sub$sampleID == "K12D4"] <-  "K12D04"
VIRGO_sub$sampleID[VIRGO_sub$sampleID == "K12D6"] <-  "K12D06"
VIRGO_sub$sampleID[VIRGO_sub$sampleID == "K12D9"] <-  "K12D09"


##For the Most part the Large Pathway Doc is cleaned Now we will subset for pathway annotation


#   In future updates will need to add Nuget Score annotation

test <- VIRGO_sub %>%  dplyr::select( sampleID, rel_sample_path_count, Annotation.x)  %>%  distinct() %>% 
  group_by(Annotation.x) %>% 
  pivot_wider(names_from = sampleID, values_from = rel_sample_path_count) %>% 
  column_to_rownames(var ="Annotation.x")  

test[is.na(test)] <- 0

test_rotate <- t(test) 

test_rotate1 <- rownames_to_column(as.data.frame(test_rotate), var ="sampleID")


test_rotate1$Subject <-  str_sub(test_rotate1$sampleID, 1, 3) 


test <- as.matrix(test)
mycolor <- viridis::viridis(paletteLength)


Heatmap(test, name = "test",
        scale = "row", 
        labRow  = VIRGO_xlabels3,
        col = mycolor,
        #split = test_rotate1$Subject
)

ggplot(VIRGO_sub, aes(x=sampleID, y=Annotation.x, fill=rel_sample_path_count))+
  geom_tile() +
  facet_grid(~Subject, scales = "free")+
  theme_classic()+
  labs(x = "Days Prior to iBV") +
  scale_x_discrete(labels=xlabels2) +
  scale_fill_viridis(discrete=FALSE)

