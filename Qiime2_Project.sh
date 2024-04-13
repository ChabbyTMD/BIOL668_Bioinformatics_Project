#!/bin/bash

# Downloading Parkinson's mouse metadata.

wget \
  -O "metadata.tsv" \
  "https://data.qiime2.org/2024.2/tutorials/pd-mice/sample_metadata.tsv"

# Exploration of the mouse metadata.

qiime metadata tabulate \
  --m-input-file metadata.tsv \
  --o-visualization metadata.qzv

#   Downloading the mouse sequence manifest tsv file.
wget \
  -O "manifest.tsv" \
  "https://data.qiime2.org/2024.2/tutorials/pd-mice/manifest"

#   Download the actual sequence archive
wget \
  -O "demultiplexed_seqs.zip" \
  "https://data.qiime2.org/2024.2/tutorials/pd-mice/demultiplexed_seqs.zip"

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

  

# After demultiplexing, which sample has the lowest sequencing depth?
# recip.460.WT.HC3.D14

# What is the median sequence length?
# 5101.5

# What is the median quality score at position 125?
# From the interactive plot, the median quality score is 38

# If you are working on this tutorial alongside someone else, why does your plot look slightly different from your neighbors? If you aren’t working alongside someone else, try running this command a few times and compare the results.

# Read trimming since full length of read is of good quality
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




# How many total features remain after denoising?

# Which sample has the highest total count of features? How many sequences did that sample have prior to DADA2 denoising?

# How many samples have fewer than 4250 total features?

# Which features are observed in at least 47 samples?

# Which sample has the fewest features? How many does it have?

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

#   Diversity analysis
qiime diversity core-metrics-phylogenetic \
  --i-table ./dada2_table.qza \
  --i-phylogeny ./tree.qza \
  --m-metadata-file ./metadata.tsv \
  --p-sampling-depth 2000 \
  --output-dir ./core-metrics-results

#   Alpha diversity
qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./core-metrics-results/faiths_pd_statistics.qzv

  qiime diversity alpha-group-significance \
 --i-alpha-diversity ./core-metrics-results/evenness_vector.qza \
 --m-metadata-file ./metadata.tsv \
 --o-visualization ./core-metrics-results/evenness_statistics.qzv

 qiime longitudinal anova \
  --m-metadata-file ./core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file ./metadata.tsv \
  --p-formula 'faith_pd ~ genotype * donor_status' \
  --o-visualization ./core-metrics-results/faiths_pd_anova.qzv

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

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column cage_id \
  --o-visualization core-metrics-results/weighted-unifrac-cage-significance_disp.qzv \
  --p-method permdisp

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