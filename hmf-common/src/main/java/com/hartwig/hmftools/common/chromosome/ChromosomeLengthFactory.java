package com.hartwig.hmftools.common.chromosome;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;

public class ChromosomeLengthFactory {

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
                        .position(samSequenceRecord.getSequenceLength())
                        .build());
            }
        }

        return results;
    }

}
