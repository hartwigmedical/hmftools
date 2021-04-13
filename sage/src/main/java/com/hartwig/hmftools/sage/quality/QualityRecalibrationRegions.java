package com.hartwig.hmftools.sage.quality;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class QualityRecalibrationRegions {

    private final IndexedFastaSequenceFile refGenome;

    public QualityRecalibrationRegions(final IndexedFastaSequenceFile refGenome) {
        this.refGenome = refGenome;
    }

    public List<GenomeRegion> regions(final int sampleSize) {
        List<GenomeRegion> result = Lists.newArrayList();

        for (final SAMSequenceRecord sequenceRecord : refGenome.getSequenceDictionary().getSequences()) {
            final String contig = sequenceRecord.getSequenceName();

            if (HumanChromosome.contains(contig) && HumanChromosome.fromString(contig).isAutosome()) {
                int start = sequenceRecord.getSequenceLength() - 1_000_000 - sampleSize;
                int end = sequenceRecord.getSequenceLength() - 1_000_001;
                result.add(GenomeRegions.create(contig, start, end));
            }
        }
        return result;
    }

}
