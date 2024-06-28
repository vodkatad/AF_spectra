MODELS=("CRC0282" "CRC1307" "CRC1502")

for M in ${MODELS[@]}; do
  echo $M
  ../local/src/new_sample_skeleton_150x.sh $M
  #../local/src/new_sample_snakes_150x.sh $M
done
