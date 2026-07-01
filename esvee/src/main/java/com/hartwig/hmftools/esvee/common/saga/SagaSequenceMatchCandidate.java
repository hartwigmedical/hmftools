package com.hartwig.hmftools.esvee.common.saga;

import java.util.List;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Cigar;

public record SagaSequenceMatchCandidate(
        SagaAlignment alignment,
        List<SagaJunctionMatchInfo> queryJunctionMatches,
        List<SagaJunctionMatchInfo> sagaJunctionMatches,
        Set<SagaSequenceMatchCandidateFilter> filters)
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
