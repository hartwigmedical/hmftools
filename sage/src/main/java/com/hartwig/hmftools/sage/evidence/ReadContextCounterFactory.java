package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.quality.UltimaModelType;
import com.hartwig.hmftools.sage.quality.UltimaQualModel;
import com.hartwig.hmftools.sage.quality.UltimaRealignedQualModels;

public class ReadContextCounterFactory
{
    private static final Set<VariantTier> HIGH_COVERAGE = EnumSet.of(VariantTier.HOTSPOT, VariantTier.PANEL);

    private final SageConfig mConfig;

    public ReadContextCounterFactory(final SageConfig config)
    {
        mConfig = config;
    }

    // TODO move elsewhere
    private static boolean isCleanHpTransition(VariantReadContext readContext, UltimaQualModel qualModel)
    {
        if(qualModel == null || !qualModel.type().equals(UltimaModelType.HOMOPOLYMER_TRANSITION))
            return false;

        byte[] coreReadBases = Arrays.subsetArray(readContext.ReadBases, readContext.CoreIndexStart, readContext.CoreIndexEnd);
        if(coreReadBases.length != readContext.RefBases.length + readContext.variant().indelLength())
            return false;

        for(int i = 0; i < coreReadBases.length; ++i)
        {
            int offset = i > readContext.leftCoreLength() ? readContext.variant().indelLength() : 0;
            if(coreReadBases[i] != readContext.RefBases[i - offset])
                return false;
        }

        return true;
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
            {
                byte[] coreBases = Arrays.subsetArray(
                        candidate.readContext().ReadBases,
                        candidate.readContext().VarIndex-1, candidate.readContext().VarIndex+1);

                UltimaQualModel qualModel = qualityCalculator.createUltimaQualModel(candidate.variant(), coreBases);

                // TODO: awkward to place this here
                if(isCleanHpTransition(candidate.readContext(), qualModel))
                {
                    candidate.readContext().setUltimaRealignedQualModels(new UltimaRealignedQualModels(candidate.readContext(),
                            qualityCalculator.ultimaQualityCalculator(), Lists.newArrayList()));
                }
                else
                {
                    candidate.readContext().setUltimaRealignedQualModels(
                            qualityCalculator.createRealignedUltimaQualModels(candidate.readContext()));
                }
            }

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
