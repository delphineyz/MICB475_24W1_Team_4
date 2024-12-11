qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path /mnt/datasets/project_2/nasa/nasa_dust_manifest.tsv \
  --output-path demux_seqs.qza

qiime demux summarize \
  --i-data demux_seqs.qza \
  --o-visualization demux_seqs.qzv

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /data/nasa_data/data_process/demux_seqs.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 262 \
  --p-trunc-len-r 233 \
  --o-representative-sequences /data/nasa_data/data_process/rep-seqs.qza \
  --o-table /data/nasa_data/data_process/table.qza \
  --o-denoising-stats /data/nasa_data/data_process/stats.qza

qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file /data/nasa_data/data_process/updated_nasa_dust_metadata.tsv \
  --o-filtered-table filtered_table.qza

qiime feature-table summarize \
  --i-table filtered_table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file /data/nasa_data/data_process/updated_nasa_dust_metadata.tsv
  
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime feature-classifier extract-reads \
 --i-sequences /mnt/datasets/silva_ref_files/silva-138-99-seqs.qza \
 --p-f-primer GTGCCAGCMGCCGCGGTAA \
 --p-r-primer GGACTACHVGGGTWTCTAAT \
 --p-trunc-len 233 \
 --o-reads ref-seqs-trimmed.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs-trimmed.qza \
  --i-reference-taxonomy /mnt/datasets/silva_ref_files/silva-138-99-tax.qza \
  --o-classifier classifier.qza

qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads /data/nasa_data/data_process/rep-seqs.qza \
  --o-classification taxonomy.qza

qiime taxa filter-table \
  --i-table filtered_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

qiime feature-table summarize \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --o-visualization dust_table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file /data/nasa_data/data_process/updated_nasa_dust_metadata.tsv

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime tools export \
   --input-path taxonomy.qza \
   --output-path /data/nasa_data/data_process/nasa_export/taxonomy_export

qiime tools export \
   --input-path rooted-tree.qza \
   --output-path /data/nasa_data/data_process/nasa_export/rooted_tree_export

qiime tools export \
   --input-path table.qza \
   --output-path /data/nasa_data/data_process/nasa_export/table_export

biom convert \
-i feature-table.biom \
--to-tsv \
-o feature-table.txt
