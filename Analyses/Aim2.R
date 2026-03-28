
# Creating the PhyloSeq Object ####

# Loading Data
library(phyloseq)
library(microbiome)
library(tidyverse)
library(ggVennDiagram)
library(indicspecies)
library(ANCOMBC)


taxonomy = read.delim('Datasets/taxonomy.tsv',
                      row.names = 1)

tree = read_tree('Datasets/tree.nwk')

counts = read.delim("Datasets/anemia-feature-table.txt", 
                    skip =1,
                    row.names =1 )

#Filter out missing and other from supplements
#Filter for 12 month olds and anemic

metadata = read.delim("Datasets/anemia_metadata.txt",
                      header = TRUE,
                      sep = "\t",
                      row.names = 1)%>%
  filter(supplement %in% c("MNP","FeSO4","None"), 
         age_months == 12,
         anemia == "anemic") %>% 
  select(supplement) %>% 
  filter(!if_any(everything(), ~ . == "Missing: Not collected"))

# Wrangling Data

taxonomy_formatted = taxonomy %>% 
  separate(col = Taxon,
           into = c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' ),
           sep=';', fill = 'right')%>%
  select(-Confidence) %>%
  as.matrix()

counts_formatted = counts %>% as.matrix()
colnames(counts_formatted) = sub("^X", "", colnames(counts_formatted))
meta_data = read.delim('Datasets/anemia_metadata.txt', row.names = 1, skip = 1, header = F)
meta_names = names(meta_data) 

# Creating the phyloseq object

ps = phyloseq(sample_data(metadata),
              otu_table(counts_formatted, taxa_are_rows = T),
              tax_table(taxonomy_formatted),
              tree)

saveRDS(ps, 'Datasets/phyloseq_taxonomy.rds')

# The Core Microbiome ####

psrare = rarefy_even_depth(ps, sample.size = 3956)
ps_rare_relab = transform(psrare, 'compositional')
ps_rare_relab_genus = tax_glom(ps_rare_relab, 'Genus')

# Subsetting the phyloseq object

supplement.MNP = subset_samples(ps_rare_relab_genus, supplement == "MNP")
supplement.FeSO4 = subset_samples(ps_rare_relab_genus, supplement == "FeSO4")
supplement.None = subset_samples(ps_rare_relab_genus, supplement == "None")

# Finding the core members

ASVs_MNP = core_members(supplement.MNP, detection = 0.001, prevalence = 0.8)
ASVs_FeSO4 = core_members(supplement.FeSO4, detection = 0.001, prevalence = 0.8)
ASVs_None = core_members(supplement.None, detection = 0.001, prevalence = 0.8)

ggVennDiagram(list(ASVs_MNP, ASVs_FeSO4, ASVs_None),
              set_size = 4,
              category.names = c("MNP", "FeSO4", "None"))

# Indicator Species Analysis ####

# Aggregating ASVs to a higher taxonomic level, converting phyloseq to relative
# abundance, and applying abundance filter

ps_phylum = tax_glom(ps, 'Phylum')
ps_relab = transform(ps_phylum, 'compositional')
ps_filt = filter_taxa(ps_relab, function(x) mean (x) > 0.001, TRUE)
otu_table = as.data.frame(otu_table(ps_filt))

set.seed(421)
indval = multipatt(t(otu_table),
                   cluster = ps_filt@sam_data$supplement,
                   control = how(nperm = 999))

summary(indval, indvalcomp = TRUE)

indval_table = as.data.frame(indval$sign)

taxonomy_f = as.data.frame(taxonomy_formatted)

signif_taxa = indval_table %>%
  left_join(taxonomy_f, by = colnames(...))

# Not sure if I am going to visualize this data. I will get back to this later. 

# Differential Abundance #### 

ps_glom = tax_glom(ps, 'Genus')

set.seed(421)
out = ancombc2(data = ps_glom,
               fix_formula = 'supplement',
               p_adj_method = 'BH',
               prv_cut = 0.5)

statistical_table = out$res

# Filter stats to include taxa that are differentially abundant between MNP and FeSO4

MNP_vs_FeSO4 = statistical_table %>%
  filter(diff_robust_supplementMNP==T) 
# Note:Not sure if this code is correct because I can't see the table. Edit if needed (See module 14)
# No supplements and control passed the stress test and have a p value less than 0.05

# Making the table

MNP_vs_FeSO4 %>%
  ggplot(aes(taxon, lfc_supplementMNP)) +
  geom_col() +
  coord_flip()


