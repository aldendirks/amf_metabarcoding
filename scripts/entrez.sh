#!/bin/bash

# Purpose:
# Script to download reference FASTA files from NCBI

# Information:
# Entrez Direct (EDirect) is a suite of utilities for interfacing with NCBI databases via UNIX command line 
# See this page for more information: https://www.ncbi.nlm.nih.gov/books/NBK179288/

# From the project root directory, usage:
#   bash download_ref.sh

####################################################################################################

# download EDirect
cd bin/ 
perl -MNet::FTP -e \
    '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1);
    $ftp->login; $ftp->binary;
    $ftp->get("/entrez/entrezdirect/edirect.tar.gz");'
gunzip -c edirect.tar.gz | tar xf -
rm edirect.tar.gz
export PATH=${PATH}:$PWD/bin/edirect
./bin/edirect/setup.sh

# number of sequences for different search parameters
# the filters "sequence from type" and "type material" are identical
# fetch Glomeromycotina type sequences
esearch -db nucleotide -query "Glomeromycotina [ORGN] AND sequence from type [FILT]" # 365
esearch -db nucleotide -query "Glomeromycotina [ORGN] AND type material [FILT]" # 365
esearch -db nucleotide -query "Glomeromycotina [ORGN] AND type material [FILT]" | \
efetch -format fasta > data/seqs/phylogenetics/input/ref_types_all.fasta

# there are 25 records that explicitly state they are holotype sequences
# there is a difference between the number that are marked as types and the number with holotype in the entry
# this suggests that 19 holotypes are not marked as such in the NCBI system and were not downloaded above
esearch -db nucleotide -query "Glomeromycotina [ORGN]" | \
efilter -query "holotype" # 25
esearch -db nucleotide -query "Glomeromycotina [ORGN] AND type material [FILT]" | \
efilter -query "holotype" # 6
esearch -db nucleotide -query "Glomeromycotina [ORGN] NOT type material [FILT]" | \
efilter -query "holotype" # 19
esearch -db nucleotide -query "Glomeromycotina [ORGN] NOT type material [FILT]" | \
efilter -query "holotype" | \
efetch -format fasta > data/seqs/phylogenetics/input/ref_undesignated_types_all.fasta

# what about paratypes?
esearch -db nucleotide -query "Glomeromycotina [ORGN]" | \
efilter -query "paratype" # 0

# what about isotypes?
# there are 8 isotypes, none of which were are designated as type material and thus were not downloaded above
esearch -db nucleotide -query "Glomeromycotina [ORGN]" | \
efilter -query "isotype" # 8
esearch -db nucleotide -query "Glomeromycotina [ORGN] NOT type material [FILT]" | \
efilter -query "isotype" # 8
esearch -db nucleotide -query "Glomeromycotina [ORGN] NOT type material [FILT]" | \
efilter -query "isotype" | \
efetch -format fasta > data/seqs/phylogenetics/input/ref_undesignated_isotypes_all.fasta

# what about epitypes?
esearch -db nucleotide -query "Glomeromycotina [ORGN]" | \
efilter -query "epitype" # 351
esearch -db nucleotide -query "Glomeromycotina [ORGN] AND type material [FILT]" | \
efilter -query "epitype" # 231
esearch -db nucleotide -query "Glomeromycotina [ORGN] NOT type material [FILT]" | \
efilter -query "epitype" # 120
esearch -db nucleotide -query "Glomeromycotina [ORGN] NOT type material [FILT]" | \
efilter -query "epitype" | \
efetch -format fasta > data/seqs/phylogenetics/input/ref_undesignated_epitypes_all.fasta

# what about extypes?

# how about searching for gene names?
# those 116 LSU sequences should be present in the ref fasta above becasue they are all from type material
esearch -db nucleotide -query "Glomeromycotina [ORGN] AND 28S rRNA [GENE]" # 2320
esearch -db nucleotide -query "Glomeromycotina [ORGN] AND type material [FILT] AND 28S rRNA [GENE]" # 116
esearch -db nucleotide -query "Glomeromycotina [ORGN] AND type material [FILT] AND 28S rRNA [GENE]" | \
efetch -format fasta > data/seqs/phylogenetics/input/ref_types_28S.fasta

# delete sequences that are shorter than 600 nucleotides
cat data/seqs/phylogenetics/input/discovered_types_lsu.fasta | \
awk '{y= i++ % 2 ; L[y]=$0; if(y==1 && length(L[1])>=600) {printf("%s\n%s\n",L[0],L[1]);}}' > \
data/seqs/phylogenetics/input/discovered_types_lsu_long.fasta

cat data/seqs/phylogenetics/input/ref_genera_lsu.fasta | \
awk '{y= i++ % 2 ; L[y]=$0; if(y==1 && length(L[1])>=600) {printf("%s\n%s\n",L[0],L[1]);}}' > \
data/seqs/phylogenetics/input/ref_genera_lsu_long.fasta

cat data/seqs/phylogenetics/input/ref_types_28S_ITSx.fasta | \
awk '{y= i++ % 2 ; L[y]=$0; if(y==1 && length(L[1])>=600) {printf("%s\n%s\n",L[0],L[1]);}}' > \
data/seqs/phylogenetics/input/ref_types_28S_ITSx_long.fasta

cat data/seqs/phylogenetics/input/ref_types_lsu.fasta | \
awk '{y= i++ % 2 ; L[y]=$0; if(y==1 && length(L[1])>=600) {printf("%s\n%s\n",L[0],L[1]);}}' > \
data/seqs/phylogenetics/input/ref_types_lsu_long.fasta

cat data/seqs/phylogenetics/input/ref_undesignated_holotypes_lsu.fasta | \
awk '{y= i++ % 2 ; L[y]=$0; if(y==1 && length(L[1])>=600) {printf("%s\n%s\n",L[0],L[1]);}}' > \
data/seqs/phylogenetics/input/ref_undesignated_holotypes_lsu_long.fasta

# download reference sequences from list of accessions
# all
esearch -db nucleotide -query "NG_060785 [accn] OR JQ048872 [accn] OR JQ048871 [accn] OR JQ048870 [accn] OR NG_060196 [accn] OR HM565946 [accn] OR HM565945 [accn] OR HM565944 [accn] OR DQ273827 [accn] OR KY555039 [accn] OR KY555040 [accn] OR KY555041 [accn] OR KJ944321 [accn] OR KJ944322 [accn] OR KJ944323 [accn] OR JN971072 [accn] OR JN971073 [accn] OR JN971074 [accn] OR JN971069 [accn] OR JN971070 [accn] OR JN971071 [accn] OR JN971076 [accn] OR JN971077 [accn] OR JN971078 [accn] OR AY541867 [accn] OR AY541859 [accn] OR AY541865 [accn] OR AY639235 [accn] OR AY639234 [accn] OR AY639233 [accn] OR JX122770 [accn] OR JX122771 [accn] OR JX122772 [accn] OR GU326339 [accn] OR GU326340 [accn] OR GU326341 [accn] OR MN384871 [accn] OR JN113035 [accn] OR JN113036 [accn] OR JN113037 [accn] OR AM183920 [accn] OR JX122774 [accn] OR JX122775 [accn] OR JX122776 [accn] OR FJ461883 [accn] OR FJ461824 [accn] OR AM158947 [accn] OR AM158948 [accn] OR FJ461842 [accn] OR KJ210826 [accn] OR KJ210827 [accn] OR KJ210828 [accn] OR HE962456 [accn] OR HE962446 [accn] OR HE962436 [accn] OR FR692348 [accn] OR FR692347 [accn] OR FM876836 [accn] OR FM876835 [accn] OR FM876834 [accn] OR FN547512 [accn] OR FN547511 [accn] OR FN547510 [accn] OR FR750371 [accn] OR FR750370 [accn] OR FR750369 [accn] OR FR750173 [accn] OR FR750172 [accn] OR FR750171 [accn] OR FR750156 [accn] OR FR750155 [accn] OR FR750154 [accn] OR FR750063 [accn] OR FM992400 [accn] OR FM992401 [accn] OR FM992402 [accn] OR HG977198 [accn] OR HG977199 [accn] OR HG977200 [accn] OR MK348926 [accn] OR MK348928 [accn] OR MK348933 [accn] OR MH560605 [accn] OR MH560606 [accn] OR MH560607 [accn] OR MK036781 [accn] OR MK036782 [accn] OR MK036783 [accn] OR MK036785 [accn] OR MK036786 [accn] OR MK036787 [accn] OR LR596344 [accn] OR LR596346 [accn] OR LR596349 [accn] OR MG459211 [accn] OR MG459212 [accn] OR MG459213 [accn] OR KF154760 [accn] OR KF154761 [accn] OR KF154762 [accn] OR KF154764 [accn] OR KF154765 [accn] OR KF154766 [accn] OR KF154768 [accn] OR KF154769 [accn] OR KF154770 [accn] OR AJ972464 [accn] OR AJ972465 [accn] OR AJ972466 [accn] OR KX345938 [accn] OR KX345939 [accn] OR KX345940 [accn] OR KY630235 [accn] OR KY630236 [accn] OR KY630237 [accn] OR KX758115 [accn] OR KX758116 [accn] OR KX758117 [accn] OR KX758121 [accn] OR KX758122 [accn] OR KX758123 [accn] OR MK903005 [accn] OR MK903006 [accn] OR MK903007 [accn] OR KJ564145 [accn] OR KJ564146 [accn] OR KJ564147 [accn] OR HG938301 [accn] OR HG938302 [accn] OR HG938303 [accn] OR KY555051 [accn] OR KY555052 [accn] OR KY555053 [accn] OR MK875630 [accn] OR MK875631 [accn] OR MK875632 [accn] OR MK570912 [accn] OR MK570913 [accn] OR MK570914 [accn] OR KF060318 [accn] OR KF060320 [accn] OR KF060321 [accn] OR KF060324 [accn] OR KF060325 [accn] OR KF060326 [accn] OR AJ849468 [accn] OR AJ849467 [accn] OR KX529097 [accn] OR KX529098 [accn] OR KX529099 [accn] OR MN306205 [accn] OR MN306206 [accn] OR MN306207 [accn] OR HG798895 [accn] OR HG798896 [accn] OR HG798897 [accn] OR KJ564166 [accn] OR KJ564167 [accn] OR KJ564168 [accn] OR KJ850181 [accn] OR KJ850182 [accn] OR KJ850183 [accn] OR KT444717 [accn] OR KT444718 [accn] OR KT444719 [accn] OR KT444708 [accn] OR KT444709 [accn] OR KT444710 [accn] OR KT444712 [accn] OR KT444713 [accn] OR KT444714 [accn] OR HG518628 [accn] OR HG518629 [accn] OR MG710517 [accn] OR MG710518 [accn] OR MG710519 [accn] OR MN384870 [accn] OR MN384872 [accn] OR MN130951 [accn] OR MN130952 [accn] OR MN130953 [accn] OR MN130955 [accn] OR MN130956 [accn] OR MN130957 [accn] OR MN130958 [accn] OR MN130959 [accn] OR MN130960 [accn] OR MH560600 [accn] OR MH560601 [accn] OR MH560602 [accn] OR MG836662 [accn] OR MG836663 [accn] OR MG836664 [accn] OR MG836659 [accn] OR MG836660 [accn] OR MG836661 [accn] OR FR865446 [accn] OR FR865447 [accn] OR FR865448 [accn] OR FR865444 [accn] OR FR865445 [accn] OR FR865449 [accn] OR FR865450 [accn] OR FR865451 [accn] OR KY630229 [accn] OR KY630230 [accn] OR KY630231 [accn] OR KJ564133 [accn] OR KJ564134 [accn] OR KJ564135 [accn] OR KJ564139 [accn] OR KJ564140 [accn] OR KJ564141 [accn] OR KJ564151 [accn] OR KJ564152 [accn] OR KJ564153 [accn] OR KJ564157 [accn] OR KJ564158 [accn] OR KJ564159 [accn] OR KJ564163 [accn] OR KJ564164 [accn] OR KJ564165 [accn] OR FR692349 [accn] OR FR692350 [accn] OR FR692351 [accn] OR JF439094 [accn] OR JF439095 [accn] OR JF439096 [accn] OR FN547535 [accn] OR FN547539 [accn] OR FN547540 [accn] OR FR750021 [accn] OR FR750022 [accn] OR FR750023 [accn] OR FR750052 [accn] OR FR750053 [accn] OR FR750054 [accn] OR FM876831 [accn] OR FM876832 [accn] OR FN547547 [accn] OR FN547548 [accn] OR FN547549 [accn] OR MN644490 [accn] OR MN644489 [accn] OR MN644490 [accn] OR KY630227 [accn] OR KY630228 [accn] OR KJ850198 [accn] OR KJ850199 [accn] OR KJ850200 [accn] OR KJ850201 [accn] OR KJ850202 [accn] OR KJ850203 [accn] OR KJ850195 [accn] OR KJ850196 [accn] OR KJ850197 [accn] OR KJ850186 [accn] OR KJ850187 [accn] OR KJ850188 [accn] OR FM865542 [accn] OR FM865544 [accn] OR FR750071 [accn] OR FR750072 [accn] OR FR750073 [accn] OR HE817871 [accn] OR FM865577 [accn] OR FR750372 [accn] OR HG969390 [accn] OR HG969391 [accn] OR HG969392 [accn] OR FR750186 [accn] OR FR750189 [accn] OR FR750190 [accn] OR HG964396 [accn] OR HG964397 [accn] OR HG964398 [accn] OR KY362436 [accn] OR KY362437 [accn] OR KY362438 [accn] OR LS974594 [accn] OR LS974595 [accn] OR LS974596 [accn]" | \
efetch -format fasta > data/seqs/phylogenetics/input/references_all.fasta

# pLSU
esearch -db nucleotide -query "NG_060785 [accn] OR JQ048872 [accn] OR JQ048871 [accn] OR JQ048870 [accn] OR NG_060196 [accn] OR HM565946 [accn] OR HM565945 [accn] OR HM565944 [accn] OR DQ273827 [accn] OR KY555039 [accn] OR KY555040 [accn] OR KY555041 [accn] OR KJ944321 [accn] OR KJ944322 [accn] OR KJ944323 [accn] OR JN971072 [accn] OR JN971073 [accn] OR JN971074 [accn] OR JN971069 [accn] OR JN971070 [accn] OR JN971071 [accn] OR JN971076 [accn] OR JN971077 [accn] OR JN971078 [accn] OR AY541867 [accn] OR AY541859 [accn] OR AY541865 [accn] OR AY639235 [accn] OR AY639234 [accn] OR AY639233 [accn] OR JX122770 [accn] OR JX122771 [accn] OR JX122772 [accn] OR GU326339 [accn] OR GU326340 [accn] OR GU326341 [accn] OR MN384871 [accn] OR JN113035 [accn] OR JN113036 [accn] OR JN113037 [accn] OR AM183920 [accn] OR JX122774 [accn] OR JX122775 [accn] OR JX122776 [accn] OR FJ461883 [accn] OR FJ461824 [accn] OR AM158947 [accn] OR AM158948 [accn] OR FJ461842 [accn] OR KJ210826 [accn] OR KJ210827 [accn] OR KJ210828 [accn]" | \
efetch -format fasta > data/seqs/phylogenetics/input/references_pLSU.fasta

# pSSU-ITS-pLSU
esearch -db nucleotide -query "HE962456 [accn] OR HE962446 [accn] OR HE962436 [accn] OR FR692348 [accn] OR FR692347 [accn] OR FM876836 [accn] OR FM876835 [accn] OR FM876834 [accn] OR FN547512 [accn] OR FN547511 [accn] OR FN547510 [accn] OR FR750371 [accn] OR FR750370 [accn] OR FR750369 [accn] OR FR750173 [accn] OR FR750172 [accn] OR FR750171 [accn] OR FR750156 [accn] OR FR750155 [accn] OR FR750154 [accn] OR FR750063 [accn] OR FM992400 [accn] OR FM992401 [accn] OR FM992402 [accn] OR HG977198 [accn] OR HG977199 [accn] OR HG977200 [accn] OR MK348926 [accn] OR MK348928 [accn] OR MK348933 [accn] OR MH560605 [accn] OR MH560606 [accn] OR MH560607 [accn] OR MK036781 [accn] OR MK036782 [accn] OR MK036783 [accn] OR MK036785 [accn] OR MK036786 [accn] OR MK036787 [accn] OR LR596344 [accn] OR LR596346 [accn] OR LR596349 [accn] OR MG459211 [accn] OR MG459212 [accn] OR MG459213 [accn] OR KF154760 [accn] OR KF154761 [accn] OR KF154762 [accn] OR KF154764 [accn] OR KF154765 [accn] OR KF154766 [accn] OR KF154768 [accn] OR KF154769 [accn] OR KF154770 [accn] OR AJ972464 [accn] OR AJ972465 [accn] OR AJ972466 [accn] OR KX345938 [accn] OR KX345939 [accn] OR KX345940 [accn] OR KY630235 [accn] OR KY630236 [accn] OR KY630237 [accn] OR KX758115 [accn] OR KX758116 [accn] OR KX758117 [accn] OR KX758121 [accn] OR KX758122 [accn] OR KX758123 [accn] OR MK903005 [accn] OR MK903006 [accn] OR MK903007 [accn] OR KJ564145 [accn] OR KJ564146 [accn] OR KJ564147 [accn] OR HG938301 [accn] OR HG938302 [accn] OR HG938303 [accn] OR KY555051 [accn] OR KY555052 [accn] OR KY555053 [accn] OR MK875630 [accn] OR MK875631 [accn] OR MK875632 [accn] OR MK570912 [accn] OR MK570913 [accn] OR MK570914 [accn] OR KF060318 [accn] OR KF060320 [accn] OR KF060321 [accn] OR KF060324 [accn] OR KF060325 [accn] OR KF060326 [accn] OR AJ849468 [accn] OR AJ849467 [accn] OR KX529097 [accn] OR KX529098 [accn] OR KX529099 [accn] OR MN306205 [accn] OR MN306206 [accn] OR MN306207 [accn] OR HG798895 [accn] OR HG798896 [accn] OR HG798897 [accn] OR KJ564166 [accn] OR KJ564167 [accn] OR KJ564168 [accn] OR KJ850181 [accn] OR KJ850182 [accn] OR KJ850183 [accn] OR KT444717 [accn] OR KT444718 [accn] OR KT444719 [accn] OR KT444708 [accn] OR KT444709 [accn] OR KT444710 [accn] OR KT444712 [accn] OR KT444713 [accn] OR KT444714 [accn] OR HG518628 [accn] OR HG518629 [accn] OR MG710517 [accn] OR MG710518 [accn] OR MG710519 [accn] OR MN384870 [accn] OR MN384872 [accn] OR MN130951 [accn] OR MN130952 [accn] OR MN130953 [accn] OR MN130955 [accn] OR MN130956 [accn] OR MN130957 [accn] OR MN130958 [accn] OR MN130959 [accn] OR MN130960 [accn] OR MH560600 [accn] OR MH560601 [accn] OR MH560602 [accn] OR MG836662 [accn] OR MG836663 [accn] OR MG836664 [accn] OR MG836659 [accn] OR MG836660 [accn] OR MG836661 [accn] OR FR865446 [accn] OR FR865447 [accn] OR FR865448 [accn] OR FR865444 [accn] OR FR865445 [accn] OR FR865449 [accn] OR FR865450 [accn] OR FR865451 [accn] OR KY630229 [accn] OR KY630230 [accn] OR KY630231 [accn] OR KJ564133 [accn] OR KJ564134 [accn] OR KJ564135 [accn] OR KJ564139 [accn] OR KJ564140 [accn] OR KJ564141 [accn] OR KJ564151 [accn] OR KJ564152 [accn] OR KJ564153 [accn] OR KJ564157 [accn] OR KJ564158 [accn] OR KJ564159 [accn] OR KJ564163 [accn] OR KJ564164 [accn] OR KJ564165 [accn] OR FR692349 [accn] OR FR692350 [accn] OR FR692351 [accn] OR JF439094 [accn] OR JF439095 [accn] OR JF439096 [accn] OR FN547535 [accn] OR FN547539 [accn] OR FN547540 [accn] OR FR750021 [accn] OR FR750022 [accn] OR FR750023 [accn] OR FR750052 [accn] OR FR750053 [accn] OR FR750054 [accn] OR FM876831 [accn] OR FM876832 [accn] OR FN547547 [accn] OR FN547548 [accn] OR FN547549 [accn] OR MN644490 [accn] OR MN644489 [accn] OR MN644490 [accn] OR KY630227 [accn] OR KY630228 [accn] OR KJ850198 [accn] OR KJ850199 [accn] OR KJ850200 [accn] OR KJ850201 [accn] OR KJ850202 [accn] OR KJ850203 [accn] OR KJ850195 [accn] OR KJ850196 [accn] OR KJ850197 [accn] OR KJ850186 [accn] OR KJ850187 [accn] OR KJ850188 [accn] OR FM865542 [accn] OR FM865544 [accn] OR FR750071 [accn] OR FR750072 [accn] OR FR750073 [accn] OR HE817871 [accn] OR FM865577 [accn] OR FR750372 [accn] OR HG969390 [accn] OR HG969391 [accn] OR HG969392 [accn] OR FR750186 [accn] OR FR750189 [accn] OR FR750190 [accn] OR HG964396 [accn] OR HG964397 [accn] OR HG964398 [accn] OR KY362436 [accn] OR KY362437 [accn] OR KY362438 [accn] OR LS974594 [accn] OR LS974595 [accn] OR LS974596 [accn]" | \
efetch -format fasta > data/seqs/phylogenetics/input/references_pSSU-ITS-pLSU.fasta

# sp. nov. material
esearch -db nucleotide -query "Glomeromycotina [ORGN] NOT type material [FILT]" | \
efilter -query "sp. nov." | \
efetch -format fasta > data/seqs/phylogenetics/input/test.fasta