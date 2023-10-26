#CRC2566LMO-L5.pass.vcf.gz  CRC2573LMO-L8.pass.vcf.gz CRC2608PRO-L6.pass.vcf.gz

table_annovar.pl CRC2566LMO-L5.pass.vcf.gz /mnt/trcanmed/snaketree/task/variant_annotations/dataset/annovar/hg38/humandb/ --otherinfo -buildver hg38 -out CRC2566LMO-L5 -remove -protocol refGene,avsnp150,cosmic87_coding -vcfinput -operation g,f,f -nastring . -polish &> CRC2566LMO-L5.multianno.log
table_annovar.pl CRC2573LMO-L8.pass.vcf.gz /mnt/trcanmed/snaketree/task/variant_annotations/dataset/annovar/hg38/humandb/ --otherinfo -buildver hg38 -out CRC2573LMO-L8 -remove -protocol refGene,avsnp150,cosmic87_coding -vcfinput -operation g,f,f -nastring . -polish &> CRC2573LMO-L8.multianno.log
table_annovar.pl CRC2608PRO-L6.pass.vcf.gz /mnt/trcanmed/snaketree/task/variant_annotations/dataset/annovar/hg38/humandb/ --otherinfo -buildver hg38 -out CRC2608PRO-L6 -remove -protocol refGene,avsnp150,cosmic87_coding -vcfinput -operation g,f,f -nastring . -polish &> CRC2608PRO-L6.multianno.log
