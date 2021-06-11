package com.hartwig.hmftools.sage.read;

import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.variant.SageVariantTier;

import org.jetbrains.annotations.NotNull;

public class ReadContextCounterFactory
{
    private static final Set<SageVariantTier> HIGH_COVERAGE = EnumSet.of(SageVariantTier.HOTSPOT, SageVariantTier.PANEL);

    private final SageConfig mConfig;
    private final Map<String, QualityRecalibrationMap> mQualityRecalibrationMap;

    public ReadContextCounterFactory(final SageConfig config, final Map<String, QualityRecalibrationMap> qualityRecalibrationMap)
    {
        mConfig = config;
        mQualityRecalibrationMap = qualityRecalibrationMap;
    }

    public List<ReadContextCounter> create(@NotNull final String sample, @NotNull final List<Candidate> candidates)
    {
        return candidates.stream()
                .map(x -> new ReadContextCounter(sample,
                        x.variant(),
                        x.readContext(),
                        mQualityRecalibrationMap.get(sample),
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
