#!/usr/bin/env bash
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AUX="${REPO_ROOT}/scripts"

REF=PPA.final.genome.nochrY.fna
QRY=Ppar.chr.rename.fa
PREFIX=PPA.Ppar

minimap2 -t 30 -I 20G -ax asm5 --eqx ${REF} ${QRY} | samtools sort -@ 8 -m 4G -o ${PREFIX}.bam
samtools index -@ 10 ${PREFIX}.bam

syri \
  -c "${PREFIX}.bam" \
  -r ${REF} \
  -q ${QRY} \
  -F B \
  --nc 32 \
  --prefix "${PREFIX}"

#filter syri results
grep -v 'SNP' ${PREFIX}syri.out > ${PREFIX}syri.noSNP.out
awk -F'\t' '$(NF-2)=="-" {print}' ${PREFIX}syri.noSNP.out > ${PREFIX}syri.noSNP.filter.out
awk '$1==$6 {print}' ${PREFIX}syri.noSNP.filter.out > ${PREFIX}syri.noSNP.filter.samechr.out
awk -F '\t' '$3 - $2 + 1 >= 100000' ${PREFIX}syri.noSNP.filter.samechr.out > ${PREFIX}syri.noSNP.filter.samechr.100k.out

#In IGV, separately examine the read coverage of Ppa_T2T and PpoHapG at the structural variant sites in ${PREFIX}syri.noSNP.filter.samechr.100k.out.

# input=ChrE1.TRANS.out
# chr=ChrE1
# eventid=TRANS125
# prefix=PPA_ChrE1_TRANS125_assembly_error

input=ChrE1.INV.out
chr=ChrE1
eventid=INV106
prefix=PPA_ChrE1_INV106_assembly_error_v2

# input=ChrB3.INV1.out
# chr=ChrB3
# eventid=INV89
# prefix=PPA_ChrB3_INV89_assembly_error

# input=ChrX.INV123.out
# chr=ChrX
# eventid=INV123
# prefix=PPA_ChrX_INV123_assembly_error

# input=ChrX.INV116.out
# chr=ChrX
# eventid=INV116
# prefix=PPA_ChrX_INV116_assembly_error

# input=ChrC2.INV96.out
# chr=ChrC2
# eventid=INV96
# prefix=PPA_ChrC2_INV96_assembly_error

python "${AUX}/plot_syri_assembly_error_v3.py" \
  --syri ${input} \
  --syri-all ${PREFIX}syri.noSNP.filter.samechr.100k.out \
  --t2t-chr ${chr} \
  --event-id ${eventid} \
  --hifi-depth GCI_hifi_10kb_avg.depth \
  --ont-depth GCI_ont_10kb_avg.depth \
  --repeats huabao_carnivora.filteredRepeats.bed \
  --sd PPA.biser.sort.merged.1.bed \
  --pur PPA.yuandongbao.PUR.bed \
  --genes PPA.gene.new.bed \
  --mask-pur-synteny \
  --pur-mask-min-len 10000 \
  --pur-mask-pad 0 \
  --pur-mask-style hatch \
  --out ${prefix} \
  --start 39500000 \
  --end 42900000
