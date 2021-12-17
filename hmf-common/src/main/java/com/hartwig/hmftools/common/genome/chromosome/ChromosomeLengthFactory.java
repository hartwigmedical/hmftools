package com.hartwig.hmftools.common.genome.chromosome;

import java.util.List;
import java.util.Map;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public final class ChromosomeLengthFactory {

    private ChromosomeLengthFactory() {
    }

    @NotNull
    public static <T extends GenomePosition> Map<Chromosome, ChromosomeLength> create(int windowSize,
            @NotNull final ListMultimap<Chromosome, T> position) {
        final Map<Chromosome, ChromosomeLength> result = Maps.newHashMap();
        for (Chromosome chromosome : position.keySet()) {
            List<T> chromosomeList = position.get(chromosome);
            if (!chromosomeList.isEmpty()) {
                T last = chromosomeList.get(chromosomeList.size() - 1);
                final int max = last.position() + windowSize - 1;
                final String contig = last.chromosome();
                result.put(chromosome, ImmutableChromosomeLength.builder().chromosome(contig).length(max).build());
            }
        }

        return result;
    }

    @NotNull
    public static List<ChromosomeLength> create(@NotNull final SAMFileHeader header) {
        return create(header.getSequenceDictionary());
    }

    @NotNull
    public static List<ChromosomeLength> create(@NotNull final SAMSequenceDictionary dictionary) {
        final List<ChromosomeLength> results = Lists.newArrayList();

        for (final SAMSequenceRecord samSequenceRecord : dictionary.getSequences()) {
            final String sequenceName = samSequenceRecord.getSequenceName();
            if (HumanChromosome.contains(sequenceName)) {
                results.add(ImmutableChromosomeLength.builder()
                        .chromosome(sequenceName)
                        .length(samSequenceRecord.getSequenceLength())
                        .build());
            }
        }

        return results;
    }
}
