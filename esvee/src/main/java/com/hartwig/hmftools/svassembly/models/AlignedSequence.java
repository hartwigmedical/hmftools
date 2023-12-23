package com.hartwig.hmftools.svassembly.models;

import java.util.List;

import org.apache.commons.lang3.tuple.Pair;

public interface AlignedSequence extends Sequence
{
    List<Alignment> getAlignmentBlocks();

    default int getMappedLength()
    {
        return getAlignmentBlocks().stream()
                .mapToInt(block -> block.isUnmapped() ? 0 : block.Length)
                .sum();
    }

    default int getMappedSectionCount()
    {
        return (int) getAlignmentBlocks().stream()
                .filter(block -> !block.isUnmapped())
                .count();
    }

    default boolean includes(final Pair<String, Integer> position)
    {
        return includes(position.getLeft(), position.getRight());
    }

    default boolean includes(final String chromosome, final int position)
    {
        return getAlignmentBlocks().stream()
                .anyMatch(block -> block.includes(chromosome, position));
    }
}
