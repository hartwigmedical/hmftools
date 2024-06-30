package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.quality.QualityCalculator;

public class ReadContextCounterFactory
{
    private static final Set<VariantTier> HIGH_COVERAGE = EnumSet.of(VariantTier.HOTSPOT, VariantTier.PANEL);

    private final SageConfig mConfig;

    public ReadContextCounterFactory(final SageConfig config)
    {
        mConfig = config;
    }

    public List<ReadContextCounter> create(
            final List<Candidate> candidates, final SageConfig config, final QualityCalculator qualityCalculator, final String sampleId)
    {
        List<ReadContextCounter> readCounters = Lists.newArrayListWithExpectedSize(candidates.size());
        boolean isReferenceSample = config.ReferenceIds.contains(sampleId);

        int readId = 0;

        for(Candidate candidate : candidates)
        {
            if(qualityCalculator.ultimaEnabled())
                candidate.readContext().setUltimaQualModel(qualityCalculator.createUltimaQualModel(candidate.variant()));

            try
            {
                readCounters.add(new ReadContextCounter(
                        readId++, candidate.readContext(), candidate.tier(), maxCoverage(candidate), candidate.minNumberOfEvents(),
                        config, qualityCalculator, sampleId, isReferenceSample));
            }
            catch(Exception e)
            {
                SG_LOGGER.error("var({}) error building counter from readContext: {}",
                        candidate.readContext().variant(), candidate.readContext());
                e.printStackTrace();
                System.exit(1);
            }
        }

        return readCounters;
    }

    private int maxCoverage(final Candidate candidate)
    {
        return HIGH_COVERAGE.contains(candidate.tier()) ? mConfig.MaxReadDepthPanel : mConfig.MaxReadDepth;
    }
}
