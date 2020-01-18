#(base) [CRC1307_clones_biod]egrassi@hactarlogin$ ln -s ../CRC1307_simul/CRC1307LMO_segments.txt CRC1307_02_0_segments.txt
s=(CRC1307_08_0 CRC1307_09_0 CRC1307_02_1_A CRC1307_02_1_B CRC1307_02_1_E CRC1307_08_1_B CRC1307_08_1_D CRC1307_08_1_E CRC1307_09_1_B CRC1307_09_1_C CRC1307_09_1_E)  

# are these segments with or without recalibration? without:
#[egrassi@occam ~]>ls -al bit/prj/snakegatk/dataset/CRC1307/sequenza/CRC1307LMO/CRC1307LMO_segments.txt 
#-rw-r--r-- 1 egrassi egrassi Group 177381 Dec 20 10:35 bit/prj/snakegatk/dataset/CRC1307/sequenza/CRC1307LMO/CRC1307LMO_segments.txt
#[egrassi@occam ~]>

for item in ${s[*]}
do
    ln -s ../CRC1307_simul/CRC1307LMO_segments.txt ${item}_segments.txt
done
