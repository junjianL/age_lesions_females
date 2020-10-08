
## Get normalized meth values

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131013/suppl/GSE131013%5Fnormalized%5Fmatrix%2Etxt%2Egz

## Get metadata

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131013/matrix/GSE131013_series_matrix.txt.gz

## Filter to get sample information

zcat GSE131013_series_matrix.txt.gz | grep Sample_description > metadata.txt
zcat GSE131013_series_matrix.txt.gz | grep Sample_characteristics_ch1 >> metadata.txt
