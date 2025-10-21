package com.hartwig.hmftools.cobalt.metrics;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.ChromosomeData;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class PartitionBuilder
{
    private final int PartitionSize;
    private final int WindowSize;

    public static ListMultimap<Chromosome, Partition> partitionChromosomes(List<ChromosomeData> chromosomes, int partitionSize, int windowSize)
    {
        ListMultimap<Chromosome, Partition> partitions = ArrayListMultimap.create();
        for (ChromosomeData chromosome : chromosomes)
        {
            PartitionBuilder builder = new PartitionBuilder(partitionSize,windowSize);
            partitions.putAll(HumanChromosome.fromString(chromosome.Name), builder.createPartitions(chromosome));
        }
        return partitions;
    }

    public PartitionBuilder(final int partitionSize, final int windowSize)
    {
        PartitionSize = partitionSize;
        WindowSize = windowSize;
    }

    public List<Partition> createPartitions(ChromosomeData chromosome)
    {
        int numberOfPartitions = chromosome.Length / PartitionSize;
        List<Partition> partitions = new ArrayList<>(numberOfPartitions);
        for (int i = 0; i < numberOfPartitions; i++)
        {
            int start = i * PartitionSize + 1;
            int end = start + PartitionSize - 1;
            partitions.add(new Partition(chromosome.Name, start, end, WindowSize));
        }
        return partitions;
    }
}
