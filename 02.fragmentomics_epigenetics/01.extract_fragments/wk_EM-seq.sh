#!/bin/bash

NUM=8
input_file=$1

Tmp_fifo=/tmp/$$.fifo
mkfifo $Tmp_fifo
exec 8<> $Tmp_fifo
rm -f $Tmp_fifo

for (( i=1; i<=NUM; i++ )); do
    echo >&8
done

while read sample group; do
    read -u 8
    {
	zcat $sample.final.vcf.gz | perl -alne 'next if /^#/; print join "_", $F[0], $F[1], $F[3], $F[4]' > ../00.mutation/$sample.txt
        input_snp="../00.mutation/$sample.txt"
        input_bam="../00.bam/$sample.filter.bam"
        output_prefix="$sample"
        minAF=0
        extend=1000

        if [[ ! -f $input_snp || ! -f $input_bam ]]; then
            echo "Error: Missing required file for $sample"
            echo >&8
            continue
        fi

        [ -s $output_prefix.sequenced.alleles ] && mv $output_prefix.sequenced.alleles $output_prefix.sequenced.alleles.bak
        index=1

        while read chr locus ref mut extra; do
            echo -en "\rLoading $index: $chr:$locus" >/dev/stderr
            let index=$index+1

            let left=$locus-$extend
            let right=$locus+$extend

            samtools view "$input_bam" "$chr:$left-$right" | \
                perl process.sam_EM.pl - "$chr:$locus" "$locus,$ref,$mut" 1
        done < "$input_snp" > "$output_prefix.sequenced.alleles"

        echo "Genotype extraction completed for $sample."
	
	mv $sample.sequenced.alleles $sample.sequenced.alleles.new1	

        perl select.pl $sample.sequenced.alleles.new1 > $sample.sequenced.alleles.new2

        perl rm_more_loci.pl $sample.sequenced.alleles.new2 > $sample.sequenced.alleles

        perl extract.info.pl $sample.sequenced.alleles $sample

        R --slave --args $sample Wt Mut < plot.size.R

        echo >&8
    } &
done < "$input_file"

wait

exec 8>&-

