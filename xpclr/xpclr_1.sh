vcftools --vcf td_filter2.vcf.recode.vcf --recode --recode-INFO-all --stdout  --keep as.txt   > td_snp_pass_as.vcf
vcftools --vcf td_filter2.vcf.recode.vcf --recode --recode-INFO-all --stdout  --keep es.txt   > td_snp_pass_es.vcf

plink --vcf td_snp_pass_es.vcf --chr chr1A --out  Chr1A.popes --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order
plink --vcf td_snp_pass_as.vcf --chr chr1A --out  Chr1A.popas --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order

plink --vcf td_snp_pass_es.vcf --chr chr2A --out  Chr2A.popes --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order
plink --vcf td_snp_pass_as.vcf --chr chr2A --out  Chr2A.popas --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order

plink --vcf td_snp_pass_es.vcf --chr chr3A --out  Chr3A.popes --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order
plink --vcf td_snp_pass_as.vcf --chr chr3A --out  Chr3A.popas --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order

plink --vcf td_snp_pass_es.vcf --chr chr4A --out  Chr4A.popes --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order
plink --vcf td_snp_pass_as.vcf --chr chr4A --out  Chr4A.popas --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order

plink --vcf td_snp_pass_es.vcf --chr chr5A --out  Chr5A.popes --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order
plink --vcf td_snp_pass_as.vcf --chr chr5A --out  Chr5A.popas --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order

plink --vcf td_snp_pass_es.vcf --chr chr6A --out  Chr6A.popes --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order
plink --vcf td_snp_pass_as.vcf --chr chr6A --out  Chr6A.popas --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order

plink --vcf td_snp_pass_es.vcf --chr chr7A --out  Chr7A.popes --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order
plink --vcf td_snp_pass_as.vcf --chr chr7A --out  Chr7A.popas --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order


plink --vcf td_snp_pass_es.vcf --chr chr1B --out  Chr1B.popes --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order
plink --vcf td_snp_pass_as.vcf --chr chr1B --out  Chr1B.popas --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order

plink --vcf td_snp_pass_es.vcf --chr chr2B --out  Chr2B.popes --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order
plink --vcf td_snp_pass_as.vcf --chr chr2B --out  Chr2B.popas --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order

plink --vcf td_snp_pass_es.vcf --chr chr3B --out  Chr3B.popes --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order
plink --vcf td_snp_pass_as.vcf --chr chr3B --out  Chr3B.popas --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order

plink --vcf td_snp_pass_es.vcf --chr chr4B --out  Chr4B.popes --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order
plink --vcf td_snp_pass_as.vcf --chr chr4B --out  Chr4B.popas --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order

plink --vcf td_snp_pass_es.vcf --chr chr5B --out  Chr5B.popes --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order
plink --vcf td_snp_pass_as.vcf --chr chr5B --out  Chr5B.popas --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order

plink --vcf td_snp_pass_es.vcf --chr chr6B --out  Chr6B.popes --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order
plink --vcf td_snp_pass_as.vcf --chr chr6B --out  Chr6B.popas --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order

plink --vcf td_snp_pass_es.vcf --chr chr7B --out  Chr7B.popes --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order
plink --vcf td_snp_pass_as.vcf --chr chr7B --out  Chr7B.popas --recode 01 transpose -output-missing-genotype 9 --allow-extra-chr --keep-allele-order


