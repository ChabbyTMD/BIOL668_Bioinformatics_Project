#!/bin/bash

# Setting up QIIME2 conda environment
# wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py38-linux-conda.yml
# conda config --set channel_priority disabled
# conda env create -n qiime2-amplicon-2024.2 --file qiime2-amplicon-2024.2-py38-linux-conda.yml
# rm qiime2-amplicon-2024.2-py38-linux-conda.yml

# Downloading Parkinson's mouse metadata.

wget \
  -O "metadata.tsv" \
  "https://data.qiime2.org/2024.2/tutorials/pd-mice/sample_metadata.tsv"

# Exploration of the mouse metadata.

qiime metadata tabulate \
  --m-input-file metadata.tsv \
  --o-visualization metadata.qzv

#   Importing data into QIIME 2

#   Downloading the mouse sequence manifest tsv file.
wget \
  -O "manifest.tsv" \
  "https://data.qiime2.org/2024.2/tutorials/pd-mice/manifest"

#   Download the actual sequence archive
wget \
  -O "demultiplexed_seqs.zip" \
  "https://data.qiime2.org/2024.2/tutorials/pd-mice/demultiplexed_seqs.zip"

unzip demultiplexed_seqs.zip
# Clean up
rm demultiplexed_seqs.zip
#   Create a QIIME artifact importing our sequences using the sequence manifest.
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path ./manifest.tsv \
  --output-path ./demux_seqs.qza

#   Visualizing demultiuplexing metrics
qiime demux summarize \
  --i-data ./demux_seqs.qza \
  --o-visualization ./demux_seqs.qzv

  
#   Question.
# After demultiplexing, which sample has the lowest sequencing depth?
# recip.460.WT.HC3.D14

# What is the median sequence length?
# 5101.5

# What is the median quality score at position 125?
# From the interactive plot, the median quality score is 38

# If you are working on this tutorial alongside someone else, why does your plot look slightly different from your neighbors? If you aren’t working alongside someone else, try running this command a few times and compare the results.

# Read trimming since full length of read is of good quality


# Sequence quality control and feature table
qiime dada2 denoise-single \
  --i-demultiplexed-seqs ./demux_seqs.qza \
  --p-trunc-len 150 \
  --o-table ./dada2_table.qza \
  --o-representative-sequences ./dada2_rep_set.qza \
  --o-denoising-stats ./dada2_stats.qza

#   Review the statistics of the denoising(trimming step)
qiime metadata tabulate \
  --m-input-file ./dada2_stats.qza  \
  --o-visualization ./dada2_stats.qzv

#   Obtain a table of the counts of each feature and each sequence.
qiime feature-table summarize \
  --i-table ./dada2_table.qza \
  --m-sample-metadata-file ./metadata.tsv \
  --o-visualization ./dada2_table.qzv



# Question

# How many total features remain after denoising?
# 287 Features remain.

# Which sample has the highest total count of features? How many sequences did that sample have prior to DADA2 denoising?
# recip.539.ASO.PD4.D14. Prior to DADA2 denoising, this sample had 5475 sequences

# How many samples have fewer than 4250 total features?
# 25 samples

# Which features are observed in at least 47 samples?
# 04c8be5a3a6ba2d70446812e99318905, ea2b0e4a93c24c6c3661cbe347f93b74 and 1ad289cd8f44e109fd95de0382c5b252 are observed in at least 47 samples.

# Which sample has the fewest features? How many does it have?
# The sample recip.460.WT.HC3.D49 has the fewest features with only 347

# Generating a phylogenetic tree for diversity analysis


#   Download reference database for the phylogentic step
wget \
  -O "sepp-refs-gg-13-8.qza" \
  "https://data.qiime2.org/2024.2/common/sepp-refs-gg-13-8.qza"


qiime fragment-insertion sepp \
	--i-representative-sequences ./dada2_rep_set.qza \
	--i-reference-database sepp-refs-gg-13-8.qza \
	--o-tree ./tree.qza \
	--o-placements ./tree_placements.qza \
	--p-threads 8

#   Alpha rarefaction
qiime diversity alpha-rarefaction \
  --i-table ./dada2_table.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./alpha_rarefaction_curves.qzv \
  --p-min-depth 10 \
  --p-max-depth 4250

# Question

# Start by opening the alpha rarefaction visualization.

#     Are all metadata columns represented in the visualization? If not, which columns were excluded and why?
#     The days_post_transplant column is not represented in the visulization because it does not contain catergorical data. Looking at the column it contains numerical data.

#     Which metric shows saturation and stabilization of the diversity?
#     observed features.

#     Which mouse genetic background has higher diversity, based on the curve? Which has shallower sampling depth?
#     From the genotype curve, the wild type mice have higher diversity. Wild type also has a shallower sampling depth.

# Now, let’s check the feature table summary.

#     What percentage of samples are lost if we set the rarefaction depth to 2500 sequences per sample?
qiime diversity alpha-rarefaction \
  --i-table ./dada2_table.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./alpha_rarefaction_curves_max2500.qzv \
  --p-min-depth 10 \
  --p-max-depth 2500



#   Diversity analysis
qiime diversity core-metrics-phylogenetic \
  --i-table ./dada2_table.qza \
  --i-phylogeny ./tree.qza \
  --m-metadata-file ./metadata.tsv \
  --p-sampling-depth 2000 \
  --output-dir ./core-metrics-results

# Question

# Where did we get the value 2000 from? Why did we pick that?
# This is the rarefaction depth that achieves saturation and preserves the highest number of samples.

#   Alpha diversity
qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./core-metrics-results/faiths_pd_statistics.qzv

  qiime diversity alpha-group-significance \
 --i-alpha-diversity ./core-metrics-results/evenness_vector.qza \
 --m-metadata-file ./metadata.tsv \
 --o-visualization ./core-metrics-results/evenness_statistics.qzv

# Question

# Is there a difference in evenness between genotype? Is there a difference in phylogenetic diversity between genotype?
# There is no difference in eveness between genotype (H=0.435, p=0.50). Additionally, there is no difference in phylogenetic diversity between genotype (H=2.347, p=0.125)

# Based on the group significance test, is there a difference in phylogenetic diversity by genotype? Is there a difference based on the donor?
# There is a marginally significant difference in phylogenetic diversity based on donor (H=3.91, p=0.047)

 qiime longitudinal anova \
  --m-metadata-file ./core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file ./metadata.tsv \
  --p-formula 'faith_pd ~ genotype * donor_status' \
  --o-visualization ./core-metrics-results/faiths_pd_anova.qzv

# Question

#     Open the unweighted UniFrac emperor plot (core-metrics-results/unweighted_unifrac_emperor.qzv) first. Can you find separation in the data? If so, can you find a metadata factor that reflects the separation? What if you used weighted UniFrac distance (core-metrics-results/weighted_unifrac_emperor.qzv)?
#     Using the unweighted UniFrac emperor plot, we can find seperation in the data while looking at the donor/donor-status and genotype and donor-status. Considering donor/donor-status the healthy and Parkinson Disease mice seperate into distinct clusters, with the exception of one outlier. While considering the genotype and  donor-status, the healthy susceptible and wildtype mice cluster into one group while the Parkinson Disease susceptible and wildtype mice cluster into their own distinct group. These findings are not different when we look at the weighted UniFrac Distance as the same patterns seem to emerge.

#     One of the major concerns in mouse studies is that sometimes differences in communities are due to natural variation in cages. Do you see clustering by cage?
# We see some clustering by cage in the weighted UniFrac distance graph. C49, C44 and C43 seem to cluster together in a particular plane while there seems to be a higher degree of dispersion for the other cage id's. The same clustering pattern with the cage id's in question is observed in the unweighted UniFrac Distance plot.



#   Beta diversity
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column donor \
  --o-visualization core-metrics-results/unweighted-unifrac-donor-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column donor \
  --o-visualization core-metrics-results/weighted-unifrac-donor-significance.qzv

# Question 

# Is there a significant effect of donor?
# There is a significant effect of donor in both the unweighted(pseudo-F=19.88, p=0.001) and weighted(pseudo-F=18.08, p=0.001) unifrac distance measures.

# From the metadata, we know that cage C31, C35, and C42 all house mice transplanted from one donor, and that cages C43, C44, and C49 are from the other. Is there a significant difference in the microbial communities between samples collected in cage C31 and C35? How about between C31 and C43? Do the results look the way you expect, based on the boxplots for donor?

# Looking at the boxplots for the weighted unifrac distance cage measures we see that there is not a significant difference in cages C31 and C35 but a not so significant difference in median microbial community composition in cages C31 and C43. In the unweighted unifrac distance we can see a significant difference in the microbial communities between C31 and C43, while there was not a significant difference in the interquartile range of C31 and C35. The results appear as what we expect in the unweighted unifrac distance plot while the weighted unifrac plot deviates from our expectations by not clearly showing the difference in microbial communities in cage C31 and C43.


qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column cage_id \
  --o-visualization core-metrics-results/unweighted-unifrac-cage-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column cage_id \
  --o-visualization core-metrics-results/weighted-unifrac-cage-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column cage_id \
  --o-visualization core-metrics-results/weighted-unifrac-cage-significance_disp.qzv \
  --p-method permdisp

# Question

# Is there a significant difference in variance for any of the cages?
# There is not a significant difference in the variance of the cages (F=1.156883,p=0.235)

qiime diversity adonis \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/unweighted_adonis.qzv \
  --p-formula genotype+donor

# Question

# If you adjust for donor in the adonis model, do you retain an effect of genotype? What percentage of the variation does genotype explain?
# In the adjusted model we retain a significant effect of genotype(p=0.011). Additionally, we see genotype explains 4% of variation in our model.

# Taxonomic classification

# Downloading pretrained classifier.
wget \
  -O "gg-13-8-99-515-806-nb-classifier.qza" \
  "https://data.qiime2.org/2024.2/common/gg-13-8-99-515-806-nb-classifier.qza"

qiime feature-classifier classify-sklearn \
  --i-reads ./dada2_rep_set.qza \
  --i-classifier ./gg-13-8-99-515-806-nb-classifier.qza \
  --o-classification ./taxonomy.qza

qiime metadata tabulate \
  --m-input-file ./taxonomy.qza \
  --o-visualization ./taxonomy.qzv      

# Create a table of representative sequences.

qiime feature-table tabulate-seqs \
  --i-data ./dada2_rep_set.qza \
  --o-visualization ./dada2_rep_set.qzv


# Question

# Find the feature, 07f183edd4e4d8aef1dcb2ab24dd7745. What is the taxonomic classification of this sequence? What’s the confidence for the assignment?
# 	k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Christensenellaceae; is the taxon of the feature with a confidence of 98.3%

# How many features are classified as g__Akkermansia?
# Only 2 features are classified as g__Akkermansia.

# Use the tabulated representative sequences to look up these features. If you blast them against NCBI, do you get the same taxonomic identifier as you obtained with q2-feature-classifier?

# Taxonomy barchart
# Visualise the different taxonomic compostion of our samples.

# Filter out samples that have a rarefaction depth lower than 2000
qiime feature-table filter-samples \
  --i-table ./dada2_table.qza \
  --p-min-frequency 2000 \
  --o-filtered-table ./table_2k.qza

# Create interactive barplot per sample 
qiime taxa barplot \
  --i-table ./table_2k.qza \
  --i-taxonomy ./taxonomy.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./taxa_barplot.qzv

# Question

# Visualize the data at level 2 (phylum level) and sort the samples by donor, then by genotype. Can you observe a consistent difference in phylum between the donors? Does this surprise you? Why or why not?
# There is a difference in the composition between the HC and the PD donors. The HC donors seem to have more phyla of batceria present than the PD donors. Interestingly samples from PD donors seem to be nearly devoid of the Verrucomicrobia and Actinobacteria phyla of bacteria. When we look at the barplots that catergorize based on genotype we notice presence of Actinobacteria bacteria in both susceptible and wildtype individuals. However we do see greater abundance of Verrucomicrobia in susceptible individuals.

# Differential abundance with ANCOM-BC

qiime feature-table filter-features \
  --i-table ./table_2k.qza \
  --p-min-frequency 50 \
  --p-min-samples 4 \
  --o-filtered-table ./table_2k_abund.qza

#   ### This Errors out due to a seemingly missing dependency libgsl.so.25
# Plugin error from composition:

#   An error was encountered while running ANCOM-BC in R (return code 1), please inspect stdout and stderr to learn more.
# # check whether there is a difference in the gut microbiome of the mice based on their donor and their genetic background.

# qiime composition ancombc \
#   --i-table ./table_2k_abund.qza \
#   --m-metadata-file ./metadata.tsv \
#   --p-formula 'donor' \
#   --o-differentials ./ancombc_donor.qza

# Commented code failed to run because the ancombc_donor artifact was not generated.

# qiime composition da-barplot \
#   --i-data ./ancombc_donor.qza \
#   --p-significance-threshold 0.001 \
#   --o-visualization da_barplot_donor.qzv

# qiime composition ancombc \
#   --i-table ./table_2k_abund.qza \
#   --m-metadata-file ./metadata.tsv \
#   --p-formula 'genotype' \
#   --o-differentials ./ancombc_genotype.qza

# qiime composition da-barplot \
#   --i-data ./ancombc_genotype.qza \
#   --p-significance-threshold 0.001 \
#   --o-visualization da_barplot_genotype.qzv

# qiime composition ancombc \
#   --i-table ./table_2k_abund.qza \
#   --m-metadata-file ./metadata.tsv \
#   --p-formula 'donor * genotype' \
#   --o-differentials ./ancombc_donor_genotype.qza

# qiime composition da-barplot \
#   --i-data ./ancombc_donor_genotype.qza \
#   --p-significance-threshold 0.001 \
#   --o-visualization da_barplot_donor_genotype.qzv


# It is not immediately clear 
#   When you open the differential abundance bar plots generated from your ANCOM-BC results, you’ll see one bar plot per column from the group(s) included in the formula parameter in the ANCOM-BC output, excluding the selected intercept(s) pulled from the reference level parameter. If a particular group is not specified in the reference level parameter, the intercept will default to the highest alphanumeric group (e.g. in alphabetical or numeric order, as applicable) within each formula term.

# Each plot visualizes features in each group compared to the intercept as LFC (log-fold change), sorted by the most relatively enriched feature to the most relatively depleted feature. As mentioned above, this visualization can be filtered by q-value (i.e., false discovery rate corrected p-value).




# Question

# Open the da-barplot visualizations for donor and genotype as the selected ANCOM-BC formula term.

#     Are there more differentially abundant features between the donors or the mouse genotype? Did you expect this result based on the beta diversity?
#  It appears that are more differentially abundant features present in the donors compared to mouse genotype.
#     Are there any features that are differentially abundant in both the donors and by genotype?
# It is not immediately clear if there are features that are differentially abundant in both the donors an by genotype

#     How do the bar plots for the combined formula (‘donor + genotype’) compare with the individual donor and mouse genotype bar plots? Are there more differentially abundant features in the individual plots or the combined?
# there appear to be more differentially abundant features in the individual plots for the donors compared to the result of the combined formula. However, the plot for the combined formula contains more differentially abundant features than the individual genotype plot.


# Taxonomic Classification
wget \
  -O "ref_seqs_v4.qza" \
  "https://data.qiime2.org/2024.2/tutorials/pd-mice/ref_seqs_v4.qza"

wget \
  -O "ref_seqs_v4.qza" \
  "https://data.qiime2.org/2024.2/tutorials/pd-mice/ref_seqs_v4.qza"

wget \
  -O "ref_tax.qza" \
  "https://data.qiime2.org/2024.2/tutorials/pd-mice/ref_tax.qza"

wget \
  -O "animal_distal_gut.qza" \
  "https://data.qiime2.org/2024.2/tutorials/pd-mice/animal_distal_gut.qza"

# Retrain the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ./ref_seqs_v4.qza \
  --i-reference-taxonomy ./ref_tax.qza \
  --i-class-weight ./animal_distal_gut.qza \
  --o-classifier ./bespoke.qza

# classify samples with our newly trained classifier.
qiime feature-classifier classify-sklearn \
  --i-reads ./dada2_rep_set.qza \
  --i-classifier ./bespoke.qza \
  --o-classification ./bespoke_taxonomy.qza

# Moving created files into new directory.
mkdir ./multi-taxonomy
mv ./taxonomy.qza ./multi-taxonomy
mv ./bespoke_taxonomy.qza ./multi-taxonomy

qiime feature-table tabulate-seqs \
  --i-data dada2_rep_set.qza \
  --i-taxonomy multi-taxonomy/ \
  --o-visualization dada2_rep_set_multi_taxonomy.qzv




# Question

# Open up the dada2_rep_set_multi_taxonomy.qzv visualization and the da_barplot_donor.qzv visualization.

#     Examine the enriched ASVs in the da_barplot_donor.qzv visualization. Are there any of these enriched ASVs that have differing taxonomic resolution in the dada2_rep_set_multi_taxonomy.qzv visualization?

# Of the enriched ASVs, the only two that differed in taxonomic resolution between the da_barplot and the dada2 visualizations were 04195686f2b70585790ec75320de0d6f and 54f7ee881a58ad84fe3f81d76968b072. The rest had similar taxonomic reslutons under both models.

#     If so, which taxonomy provided better resolution?
# The bespoke taxonomy provided better taxonomic resolution for more features overall. However the stock model had more detalied taxonomic resolution for a few samples. One of them being 54f7ee881a58ad84fe3f81d76968b072

#     Is this what we expect, based on what we learned about taxonomic classification, accuracy, and re-training earlier in the tutorial?
# this fits our expectations. Since the stock model is a good starting point that captures a majority of the taxonomic clasification for most features, having a bespoke model for our samples allows us to achive finer taxonomic resolution for some of our features.


# PCoA-based analyses

# Question
    # Open the unweighted UniFrac emperor plot and color the samples by mouse id. Click on the “animations” tab and animate using the day_post_transplant as your gradient and mouse_id as your trajectory. Do you observe any clear temporal trends based on the PCoA?
    # The trajectyories seem to be centered within the clustered, aside from one that traverses into another cluster.

    # Can we visualize change over time without an animation? What happens if you color the plot by day_post_transplant? Do you see a difference based on the day? Hint: Try changing the colormap to a sequential colormap like viridis.

# Using the controls, look at variation in cage along PCs 1, 2, and 3. What kind of patterns do you see with time along each axis?


qiime longitudinal volatility \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-file ./core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --p-state-column days_post_transplant \
  --p-individual-id-column mouse_id \
  --p-default-group-column 'donor_status' \
  --p-default-metric 'Axis 2' \
  --o-visualization ./pc_vol.qzv

# Distance-based analysis
qiime longitudinal first-distances \
  --i-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata.tsv \
  --p-state-column days_post_transplant \
  --p-individual-id-column mouse_id \
  --p-baseline 7 \
  --o-first-distances ./from_first_unifrac.qza

  qiime longitudinal volatility \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-file ./from_first_unifrac.qza \
  --p-state-column days_post_transplant \
  --p-individual-id-column mouse_id \
  --p-default-metric Distance \
  --p-default-group-column 'donor_status' \
  --o-visualization ./from_first_unifrac_vol.qzv
  
# Question
  # Based on the volatility plot, does one donor change more over time than the other? What about by genotype? Cage?
  # Ther seems to be an overall upward trend with healthy donors and a downward trend for PD donors. There also seems to be an upward trend for all cage ID's except C44 and C49 which show a downward trend. Wildtype trends upward while susceptible trends upwards for about 5 days and plateaus.


  qiime longitudinal linear-mixed-effects \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-file ./from_first_unifrac.qza \
  --p-metric Distance \
  --p-state-column days_post_transplant \
  --p-individual-id-column mouse_id \
  --p-group-columns genotype,donor \
  --o-visualization ./from_first_unifrac_lme.qzv



# Question

#     Is there a significant association between the genotype and temporal change?
# There is a significant association between genotype and temporal change(p=0.044)
#     Which genotype is more stable (has lower variation)?
# The PD genotype has lower variation
#     Is there a temporal change associated with the donor? Did you expect or not expect this based on the volatility plot results?
# There is a significant assocaition between temporal change and donor(p=0.018). This agrees with the volatility plot results that showed an upward trend with increasing days post transplant.
#     Can you find an interaction between the donor and genotype?
# There seems to be a significant association between donor and genotype (p=0.003)

# Machine-learning classifiers for predicting sample characteristics

qiime sample-classifier classify-samples \
  --i-table ./dada2_table.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column genotype_and_donor_status \
  --p-random-state 666 \
  --p-n-jobs 1 \
  --output-dir ./sample-classifier-results/



# Question

# How did we do? Just for fun, try predicting some of the other metadata columns to see how easily cage_id and other columns can be predicted.

qiime sample-classifier heatmap \
  --i-table ./dada2_table.qza \
  --i-importance ./sample-classifier-results/feature_importance.qza \
  --m-sample-metadata-file ./metadata.tsv \
  --m-sample-metadata-column genotype_and_donor_status \
  --p-group-samples \
  --p-feature-count 100 \
  --o-heatmap ./sample-classifier-results/heatmap_100-features.qzv \
  --o-filtered-table ./sample-classifier-results/filtered-table_100-features.qza


