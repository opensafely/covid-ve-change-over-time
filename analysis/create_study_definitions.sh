# Run this script to create a study_definition for each k

for i in {1..6}; do
sed -e "s;%placeholder_k%;$i;g" ./analysis/study_definition_k.py > ./analysis/study_definition_$i.py;
done;
