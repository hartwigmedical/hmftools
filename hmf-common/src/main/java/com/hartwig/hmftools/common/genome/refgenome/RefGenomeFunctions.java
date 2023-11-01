package com.hartwig.hmftools.common.genome.refgenome;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.CHR_PREFIX;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SamReader;

public final class RefGenomeFunctions
{
    public static final Logger LOGGER = LogManager.getLogger(RefGenomeFunctions.class);

    private RefGenomeFunctions() {}

    public static boolean samReaderUsesChrInContigs(final SamReader samReader)
    {
        return samReader.getFileHeader()
                .getSequenceDictionary()
                .getSequences()
                .stream()
                .anyMatch(x -> x.getSequenceName().contains(CHR_PREFIX));
    }

    public static String stripChrPrefix(final String chromosome)
    {
        if(chromosome.startsWith(CHR_PREFIX))
        {
            return chromosome.substring(CHR_PREFIX.length());
        }

        return chromosome;
    }

    public static String enforceChrPrefix(final String chromosome)
    {
        if(!chromosome.startsWith(CHR_PREFIX))
        {
            return CHR_PREFIX + chromosome;
        }

        return chromosome;
    }
}
