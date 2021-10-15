#################33
### Run Emmax
### Example code
####################

### Centroidsize Example

# Same_Chr$num.vcf.gz 
### vcf file with the CNV information for the Same direction model

# Bantu_phenotypes.ped
### ped file with phenotype information

# chr$num.kinship.out 
### kinship matrix 


for num in {1..22} # loop through the chromosomes
do
./epacts single \
 --vcf Same_Chr$num.vcf.gz --ped  Bantu_phenotypes.ped \
 --min-maf 0.0 --kin chr$num\kinship.out \
 --pheno CS --cov AGE --cov SEX --test q.emmax \
 --out Same_CS_Chr$num --run 1
done
