conda activate qiime2-2020.11
qiime tools export --input-path dada2_table.qza --output-path exported_table/
biom convert -i exported_table/feature-table.biom -o exported_table/feature-table.biom.tsv --to-tsv
conda deactivate
