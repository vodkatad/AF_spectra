#(base) [CRC1307_clones_biod]egrassi@hactarlogin$ ln -s ../CRC1307_simul/CRC1307LMO_segments.txt CRC1307_02_0_segments.txt
s=(CRC1307-02-0 CRC1307-08-0 CRC1307-09-0 CRC1307-02-1-A CRC1307-02-1-B CRC1307-02-1-E CRC1307-08-1-B CRC1307-08-1-D CRC1307-08-1-E CRC1307-09-1-B CRC1307-09-1-C CRC1307-09-1-E CRC1307-08-MA-A CRC1307-08-MA-C CRC1307-08-MA-F CRC1307-08-MC-D CRC1307-08-MC-E CRC1307-08-MC-F CRC1307-08-MI-A CRC1307-08-MI-B CRC1307-08-MI-F)  

# are these segments with or without recalibration? without:
#[egrassi@occam ~]>ls -al bit/prj/snakegatk/dataset/CRC1307/sequenza/CRC1307LMO/CRC1307LMO_segments.txt 
#-rw-r--r-- 1 egrassi egrassi Group 177381 Dec 20 10:35 bit/prj/snakegatk/dataset/CRC1307/sequenza/CRC1307LMO/CRC1307LMO_segments.txt
#[egrassi@occam ~]>

for item in ${s[*]}
do
    ln -s ../CRC1307_simul/CRC1307LMO_segments.txt ${item}.segments.txt
done
