package com.hartwig.hmftools.svprep.reads;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

public class ReadGroupState
{
    public final String ReadId;
    public final List<ExpectedRead> Reads;
    public final int PartitionCount;

    // public final Set<String> RemoteChromosomePartitions;
    // public final Set<String> ProcessedChromosomePartitions;

    public ReadGroupState(final String readId, final List<ExpectedRead> reads, final int partitionCount)
    {
        ReadId = readId;
        PartitionCount = partitionCount;
        Reads = reads;

        // RemoteChromosomePartitions = Sets.newHashSet();
        // ProcessedChromosomePartitions = Sets.newHashSet();
    }

    public static ReadGroupState formGroup(final ReadGroup readGroup)
    {
        List<ExpectedRead> reads = Lists.newArrayList();

        for(ReadRecord read : readGroup.reads())
        {
            reads.add(new ExpectedRead(
                    read.Chromosome, read.start(), read.isFirstOfPair(), false, true));

            if(read.hasSuppAlignment())
            {
                SupplementaryReadData suppData = read.supplementaryAlignment();

                reads.add(new ExpectedRead(
                        suppData.Chromosome, suppData.Position, read.isFirstOfPair(), true, false));
            }

            if(HumanChromosome.contains(read.MateChromosome) && !readGroup.hasReadMate(read))
            {
                reads.add(new ExpectedRead(
                        read.MateChromosome, read.MatePosStart, !read.isFirstOfPair(), false, false));
            }
        }

        reads.forEach(x -> x.setExpectedMatchCount(readGroup.partitionCount() - 1));
        return new ReadGroupState(readGroup.id(), reads, readGroup.partitionCount());
    }

    /*
    public ExpectedRead findMatchingRead(final ExpectedRead read)
    {
        return Reads.stream().filter(x -> x.matches(read)).findFirst().orElse(null);
    }
    */

    public String toString()
    {
        return format("id(%s) reads(%d)", ReadId, Reads.size());
    }
}
