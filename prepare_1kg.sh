#!/usr/bin/env bash
#
# This script prepares the 1000 Reference Genomes data
# by editing the FIDs.
#
# Author: Nikolas Baya (2020-09-19)

original=pop_4pop_mix_SEQ
out=ref_panel_1kg_4pop

sed -e 's/mis_pop_//g' \
    -e 's/_SEQ//g' \
    -e 's/afri_afr/AFR/g' \
    -e 's/amer_amr/AMR/g' \
    -e 's/asia_asn/ASN/g' \
    -e 's/euro_eur/EUR/g' \
    ${original}.fam \
    | sed 's/*/\t/g' \
    | awk '{ print $1,$3,$4,$5,$6,$7 }' > ${out}.fam

for suffix in {bed,bim}; do
  cp ${original}.${suffix} ${out}.${suffix}
done
