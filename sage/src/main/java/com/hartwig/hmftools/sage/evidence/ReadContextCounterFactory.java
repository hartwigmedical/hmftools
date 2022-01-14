package com.hartwig.hmftools.sage.evidence;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.common.VariantTier;

import org.jetbrains.annotations.NotNull;

public class ReadContextCounterFactory
{
    private static final Set<VariantTier> HIGH_COVERAGE = EnumSet.of(VariantTier.HOTSPOT, VariantTier.PANEL);

    private final SageConfig mConfig;

    public ReadContextCounterFactory(final SageConfig config)
    {
        mConfig = config;
    }

    public List<ReadContextCounter> create(@NotNull final String sample, @NotNull final List<Candidate> candidates)
    {
        return candidates.stream()
                .map(x -> new ReadContextCounter(sample,
                        x.variant(),
                        x.readContext(),
                        x.tier(),
                        maxCoverage(x),
                        x.minNumberOfEvents(),
                        mConfig.maxSkippedReferenceRegions(),
                        x.maxReadDepth() < mConfig.MaxRealignmentDepth))
                .collect(Collectors.toList());
    }

    private int maxCoverage(@NotNull final Candidate candidate)
    {
        return HIGH_COVERAGE.contains(candidate.tier()) || MitochondrialChromosome.contains(candidate.chromosome())
                ? mConfig.MaxReadDepthPanel : mConfig.MaxReadDepth;
    }
}
