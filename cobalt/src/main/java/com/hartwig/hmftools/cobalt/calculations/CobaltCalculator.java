package com.hartwig.hmftools.cobalt.calculations;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.CobaltConfig;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.cobalt.targeted.TargetRegions;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class CobaltCalculator
{
    private final ListMultimap<Chromosome, DepthReading> mTumorReadDepths;
    private final GenomeFilter mWindowStatuses;
    private final TargetRegions mEnricher;
    private final CobaltConfig mConfig;

    public CobaltCalculator(final ListMultimap<Chromosome, DepthReading> mTumorReadDepths, CobaltConfig config)
    {
        this.mTumorReadDepths = mTumorReadDepths;
        mWindowStatuses = new WindowStatuses(config.gcProfileData(), config.mExcludedRegions);
        mConfig = config;
        mEnricher = mConfig.targetRegionEnricher();
    }

    public ListMultimap<Chromosome, CobaltRatio> doCalculation()
    {
        CobaltCalculation tumorCalculation = new CobaltCalculation(mWindowStatuses, mConfig.RefGenVersion,mEnricher);
        mTumorReadDepths.forEach((tumorCalculation::addReading));
        return tumorCalculation.calculateRatios();
    }
}
