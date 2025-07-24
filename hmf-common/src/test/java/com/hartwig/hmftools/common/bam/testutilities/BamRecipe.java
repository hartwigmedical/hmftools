package com.hartwig.hmftools.common.bam.testutilities;

import java.util.HashMap;
import java.util.Map;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.bam.FastBamWriter;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMSequenceRecord;

public class BamRecipe
{
    private final Map<Integer,ChromosomeRegionDepths> chromosomeDepths = new HashMap<>();

    public void add(ChromosomeRegionDepths depths)
    {
        Preconditions.checkArgument(!chromosomeDepths.containsKey(depths.mChromosome), "Duplicate chromosome index");
        chromosomeDepths.put(depths.mChromosome, depths);
    }

    public void writeToBam(String outputFileName, RefGenomeSource refGenomeSource)
    {
        var header = new SAMFileHeader();
        chromosomeDepths.keySet().forEach(chromosomeIndex -> {
            String chrName = "chr" + (chromosomeIndex + 1);
            int chrLength = refGenomeSource.chromosomeLengths().get(chrName);
            header.getSequenceDictionary().addSequence(new SAMSequenceRecord(chrName, chrLength));
        });
        SAMFileWriter bamWriter = new FastBamWriter(header, outputFileName);

        chromosomeDepths.values().forEach(depths -> depths.writeToBam(bamWriter, refGenomeSource));
        bamWriter.close();
    }
}