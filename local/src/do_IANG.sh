
MODELS=('CRCUECHPRO' 'CRC2826PRO' 'CRC3023PRO')
#MODELS=("CRC0441")
for M in ${MODELS[@]}; do
    echo $M
    #../local/src/new_sample_skeleton_IANG.sh $M
    ../local/src/new_sample_snakes_V2.sh $M
done

#exit 0
