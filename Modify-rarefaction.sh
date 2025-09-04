# input is rarefaction file generated with mothur saved as a .txt file, called rarefaction.txt
# code drops -lci and -hci columns from rarefaction file for plotting in R and saves as new file called rarefaction_means.txt

awk -F'\t' 'NR==1 {
    for(i=1;i<=NF;i++) if($i ~ /^numsampled$/ || $i ~ /^0\.03-/) keep[i]
}
{
    first=1
    for(i=1;i<=NF;i++) if(i in keep) {
        if(!first) printf "\t"
        printf "%s",$i
        first=0
    }
    printf "\n"
}' rarefaction.txt > rarefaction_means.txt


