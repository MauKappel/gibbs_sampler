cp -R data non-binders
cd non-binders
for file in */*
do
   echo "$file"
   awk '{ if ($2 > 0.426) { print $1 " " $2 } }' "../data/$file" > "$file"
done