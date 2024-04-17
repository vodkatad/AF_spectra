input=ovlen2
while IFS= read -r line
do
  cat $line | bawk '$2==3'
done < "$input" > ovlen_cn3_nomsi