#!/usr/bin/env bash

a=$#
if [ $# -lt 5 ]
        then 
        printf "\n\n"   
        echo "Promoter prediction in Mycoplasmas"
        printf "\n"
        echo "Usage $0 <Output name> <Fasta sequence> <Pribnow box PSSM> <-35 box PSSM> <Stacking energies table>"
        printf "\n\n"
        exit 2  
fi

#output=`echo $1 | awk '{ sub(/\.txt/, ""); print }'`
#output2="/tmp/${output}"
#outputFolder="/tmp/${output}_results"

output=`echo $1 | awk '{ sub(/\.txt/, ""); print }'`
outputFolder="${output}_results"

mkdir $outputFolder

echo "Getting all parameters.."

./promoterNewVero $2 $5 $output $3 $4

echo "Applying Random Forests.."

Rscript ./predict_promoters.R $output

resultP="${output}_final_plus.txt"
resultM="${output}_final_minus.txt"

mv $resultP $outputFolder/
mv $resultM $outputFolder/

rm1="${output}_minus.txt"
rm2="${output}_plus.txt"

#rm $rm1
#rm $rm2

dataP="${outputFolder}/${resultP}"
dataM="${outputFolder}/${resultM}"

dataFINAL="${outputFolder}/${output}_filtered.txt"

Rscript ./retrievePromoters.R $dataP $dataM $dataFINAL
