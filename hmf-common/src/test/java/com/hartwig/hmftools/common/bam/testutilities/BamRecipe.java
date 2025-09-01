package com.hartwig.hmftools.common.bam.testutilities;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.bam.FastBamWriter;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMSequenceRecord;

public class BamRecipe
{
    private final Map<HumanChromosome, List<ChromosomeRegionDepths>> chromosomeDepths = new HashMap<>();
    private final ChromosomeLengths chromosomeLengths;

    public BamRecipe(final ChromosomeLengths chromosomeLengths)
    {
        this.chromosomeLengths = chromosomeLengths;
    }

    public void add(ChromosomeRegionDepths depths)
    {
        List<ChromosomeRegionDepths> depthsForChromosome = chromosomeDepths.getOrDefault(depths.mChromosome, new ArrayList<>());
        depthsForChromosome.add(depths);
        chromosomeDepths.put(depths.mChromosome, depthsForChromosome);
    }

    public void writeToBam(String outputFileName)
    {
        var header = new SAMFileHeader();
        chromosomeDepths.keySet().forEach(chromosomeIndex ->
        {
            String chrName = RefGenomeVersion.V38.versionedChromosome(chromosomeIndex);
            int chrLength = chromosomeLengths.chromosomeLength(chromosomeIndex);
            header.getSequenceDictionary().addSequence(new SAMSequenceRecord(chrName, chrLength));
        });
        SAMFileWriter bamWriter = new FastBamWriter(header, outputFileName);

        chromosomeDepths.values()
                .forEach(depthsList ->
                        depthsList.forEach( chromosomeRegionDepths ->
                                chromosomeRegionDepths.writeToBam(bamWriter)));
        bamWriter.close();
    }
}