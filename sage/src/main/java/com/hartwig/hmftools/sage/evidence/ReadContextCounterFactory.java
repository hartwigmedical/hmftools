package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConfig.AppendMode;
import static com.hartwig.hmftools.sage.SageConfig.isUltima;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MAX_GERMLINE_VAF_PANEL_INDEL_REPEAT_THRESHOLD_FACTOR;
import static com.hartwig.hmftools.sage.seqtech.UltimaQualModelBuilder.setReadContextUltimaModels;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.sage.quality.QualityCalculator;

import org.apache.logging.log4j.Level;

public class ReadContextCounterFactory
{
    private static final Set<VariantTier> HIGH_COVERAGE = EnumSet.of(VariantTier.HOTSPOT, VariantTier.PANEL);

    private final SageConfig mConfig;

    public ReadContextCounterFactory(final SageConfig config)
    {
        mConfig = config;
    }

    public List<ReadContextCounter> create(
            final List<Candidate> candidates, final SageConfig config, final RefGenomeInterface refGenome,
            final QualityCalculator qualityCalculator, final String sampleId)
    {
        List<ReadContextCounter> readCounters = Lists.newArrayListWithExpectedSize(candidates.size());
        boolean isReferenceSample = config.ReferenceIds.contains(sampleId);

        int readId = 0;

        for(Candidate candidate : candidates)
        {
            try
            {
                ReadContextCounter readContextCounter = new ReadContextCounter(
                        readId++, candidate.readContext(), candidate.tier(), maxCoverage(candidate), candidate.minNumberOfEvents(),
                        config, qualityCalculator, sampleId, isReferenceSample);

                if(isUltima())
                {
                    setReadContextUltimaModels(refGenome, readContextCounter);
                }

                readCounters.add(readContextCounter);
            }
            catch(Exception e)
            {
                boolean isNonStandardChromosome = !HumanChromosome.contains(candidate.variant().chromosome());
                boolean isLowBoundaryPosition = candidate.variant().position() < DEFAULT_FLANK_LENGTH * 3;

                boolean allowInvalidContexts = AppendMode || isNonStandardChromosome || isLowBoundaryPosition;

                if(!isNonStandardChromosome)
                {
                    Level logLevel = allowInvalidContexts ? Level.WARN : Level.ERROR;

                    SG_LOGGER.log(logLevel, "var({}) error building counter from readContext: {}",
                            candidate.readContext().variant(), candidate.readContext());
                }

                if(allowInvalidContexts)
                {
                    candidate.readContext().markedInvalid();

                    ReadContextCounter readContextCounter = new ReadContextCounter(
                            readId++, candidate.readContext(), candidate.tier(), maxCoverage(candidate), candidate.minNumberOfEvents(),
                            config, sampleId);

                    readCounters.add(readContextCounter);
                }
                else
                {
                    e.printStackTrace();
                    System.exit(1);
                }
            }
        }

        return readCounters;
    }

    private int maxCoverage(final Candidate candidate)
    {
        return HIGH_COVERAGE.contains(candidate.tier()) ? mConfig.MaxReadDepthPanel : mConfig.MaxReadDepth;
    }
}
