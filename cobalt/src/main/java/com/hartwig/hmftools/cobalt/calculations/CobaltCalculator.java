package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import java.util.List;

import com.google.common.base.Preconditions;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.CobaltConfig;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.cobalt.targeted.CobaltScope;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.GcMedianReadDepth;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public class CobaltCalculator
{
    private final ListMultimap<Chromosome, CobaltRatio> mRatios;
    private final List<MedianRatio> mMedianRatios;
    private final GcMedianReadDepth mTumorStats;
    private final GcMedianReadDepth mReferenceStatistics;

    public CobaltCalculator(
            final ListMultimap<HumanChromosome, DepthReading> tumourDepthReadings,
            final ListMultimap<HumanChromosome, DepthReading> referenceDepthReadings,
            CobaltConfig config)
    {
        Preconditions.checkArgument(!tumourDepthReadings.isEmpty() || !referenceDepthReadings.isEmpty());
        WindowStatuses mWindowStatuses = new WindowStatuses(config.gcProfileData(), config.excludedRegions(), config.diploidRegions());
        CobaltScope scope = config.scope();

        TumorCalculation tumorCalc = new TumorCalculation(mWindowStatuses, scope, config.genome());
        tumourDepthReadings.forEach(tumorCalc::addReading);
        ListMultimap<Chromosome, BamRatio> tumorResults = tumorCalc.calculateRatios();
        mTumorStats = tumorCalc.medianReadDepths();
        CB_LOGGER.info("Tumor sample median: {}, mean: {}", mTumorStats.medianReadDepth(), mTumorStats.meanReadDepth());

        ReferenceCalculation referenceCalc = new ReferenceCalculation(mWindowStatuses, scope, config.genome(), tumorCalc.consolidator());
        referenceDepthReadings.forEach(referenceCalc::addReading);
        ListMultimap<Chromosome, BamRatio> referenceResults = referenceCalc.calculateRatios();
        mMedianRatios = referenceCalc.medianRatios();
        mReferenceStatistics = referenceCalc.medianReadDepths();
        CB_LOGGER.info("Reference sample median: {}, mean: {}", mReferenceStatistics.medianReadDepth(), mReferenceStatistics.meanReadDepth());

        ResultsCollator collator = new ResultsCollator(config.genome());
        mRatios = collator.collateResults(tumorResults, referenceResults);
    }

    public ListMultimap<Chromosome, CobaltRatio> getCalculatedRatios()
    {
        return mRatios;
    }

    public List<MedianRatio> medianRatios()
    {
        return mMedianRatios;
    }

    public GcMedianReadDepth tumorMedianReadDepth()
    {
        return mTumorStats;
    }

    public GcMedianReadDepth referenceMedianReadDepth()
    {
        return mReferenceStatistics;
    }
}
