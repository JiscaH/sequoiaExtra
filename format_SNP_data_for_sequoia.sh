#!/bin/bash

##  recode to 1 column per SNP, coded as 0/1/2 copies of minor allele, missing = NA
if test -f $1.bed; then
  plink --bfile $1 --recode A --out TmpFileName
else
  plink --file $1 --recode A --out TmpFileName
fi

## drop columns 2--6
## drop header row
## replace missing value code NA by -9
cat TmpFileName.raw | tr -s ' ' | cut -d ' ' -f2,7- > $1.txt 
sed -i '1d;s/NA/-9/g' $1.txt  
rm TmpFileName.*


## drop rows with IDs with invalid characters
## TODO
## see e.g https://opensource.com/article/18/5/you-dont-know-bash-intro-bash-arrays
# # cat {$1}.raw | wc -l > nInd

# # for n in {1..$nInd}
# # do
  # # if [[ $ID =~ ^[A-Za-z0-9._]+$ ]]; then
      # # echo "OK"
  # # else
      # # echo "Not OK"
  # # fi
# # done
