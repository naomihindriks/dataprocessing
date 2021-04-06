#!/usr/bin/env sh
input1=$(realpath $1)
output1=$(realpath $2)
input2=$(realpath $3)
output2=$(realpath $4)

qiime tools export --input-path ${input1} --output-path exported-feature-table
biom convert -i exported-feature-table/feature-table.biom -o ${output1} --to-tsv
rm -drf exported-feature-table

qiime tools export --input-path ${input2} --output-path ${output2}