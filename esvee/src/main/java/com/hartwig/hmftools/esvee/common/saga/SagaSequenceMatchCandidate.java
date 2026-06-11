package com.hartwig.hmftools.esvee.common.saga;

import java.util.List;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Cigar;

public record SagaSequenceMatchCandidate(
        SagaAlignment alignment,
        List<List<Integer>> esveeJunctionOverlaps,
        List<List<Integer>> sagaJunctionOverlaps,
        List<String> filters)
{
    public SagaAssembly sagaAssembly()
    {
        return alignment.sagaAssembly();
    }

    public SagaVariant variant()
    {
        return sagaAssembly().variant();
    }

    public Cigar cigar()
    {
        return alignment.cigar();
    }

    public int alignScore()
    {
        return alignment.alignScore();
    }

    public boolean isAccepted()
    {
        return filters.isEmpty();
    }

    @NotNull
    @Override
    public String toString()
    {
        return String.format("SagaSequenceMatchCandidate(variant=\"%s\", cigar=%s, alignScore=%d, filters=%s)", variant(), cigar(), alignScore(), filters);
    }
}
