MODELS=("CRC0282" "CRC0327" "CRC1078" "CRC1307" "CRC1502" "CRC1599LM" "CRC1599PR" "CRC0441")
for M in ${MODELS[@]}; do
  echo -en "$M\t"
  cat ${M}/platypus_nobin_00/*prom*| cut -f 1 |sort | uniq | wc -l
done

