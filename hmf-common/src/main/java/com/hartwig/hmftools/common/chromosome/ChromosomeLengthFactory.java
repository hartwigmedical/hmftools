package com.hartwig.hmftools.common.chromosome;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;

public final class ChromosomeLengthFactory {

    @NotNull
    public static <T extends GenomePosition> Map<String, ChromosomeLength> create(int windowSize,
            @NotNull final Multimap<String, T> position) {
        final Map<String, ChromosomeLength> result = Maps.newHashMap();
        for (String chromosome : position.keySet()) {
            long max = position.get(chromosome).stream().mapToLong(x -> x.position() + windowSize - 1).max().orElse(0L);
            result.put(chromosome, ImmutableChromosomeLength.builder().chromosome(chromosome).length(max).build());
        }

        return result;
    }

    @NotNull
    public static List<ChromosomeLength> create(@NotNull final SAMFileHeader header) {
        return create(header, Sets.newHashSet("MT"));
    }

    @NotNull
    public static List<ChromosomeLength> create(@NotNull final SAMFileHeader header, @NotNull final Set<String> excludedSequences) {
        final List<ChromosomeLength> results = Lists.newArrayList();

        for (final SAMSequenceRecord samSequenceRecord : header.getSequenceDictionary().getSequences()) {
            final String chromosome = samSequenceRecord.getSequenceName();
            if (!excludedSequences.contains(chromosome)) {
                results.add(ImmutableChromosomeLength.builder()
                        .chromosome(chromosome)
                        .length(samSequenceRecord.getSequenceLength())
                        .build());
            }
        }

        return results;
    }
}
