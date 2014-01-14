

1:265,758,994-265,759,821

 -- 265,762,357

1:265,770,984-265,771,811

 -- 265,774,452

1:265,779,667-265,780,494


1:265,798,717-265,799,544


## Test Sets

    REF=~/Projects/f2-het-pilot/data/aln/Zea_mays.AGPv3.20.dna.combined.fa
    samtools view -h  ../TIL05.sorted.bam 1:265770984-265774452 | \
	  samtools view -Sbh - | \
	  samtools calmd -b - $REF > testset-1.bam


