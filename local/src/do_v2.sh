
MODELS=("CRC1599PR" "CRC0441")
#MODELS=("CRC0441")
for M in ${MODELS[@]}; do
    echo $M
    #../local/src/new_sample_skeleton_V2.sh $M
    ../local/src/new_sample_snakes_V2.sh $M
done

#exit 0

MODELS=("CRC0282" "CRC0327" "CRC1078" "CRC1307" "CRC1502" "CRC1599LM" "CRC0327")
#MODELS=("CRC0327")
#MODELS=("CRC1599LM")
for M in ${MODELS[@]}; do
  echo $M
  #../local/src/new_sample_skeleton_V2.sh $M
  ../local/src/new_sample_snakes_V2.sh $M
  #../local/src/new_sample_skeleton_V2.sh ${M}_clones_all
  #../local/src/new_sample_snakes_V2.sh ${M}_clones_all
done
