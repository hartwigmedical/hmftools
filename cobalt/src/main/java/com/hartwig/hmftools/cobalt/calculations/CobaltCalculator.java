package com.hartwig.hmftools.cobalt.calculations;

import java.io.IOException;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.CobaltConfig;
import com.hartwig.hmftools.cobalt.count.ReadDepth;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class CobaltCalculator
{
    private final ListMultimap<Chromosome, ReadDepth> mTumorReadDepths;
    private final WindowStatuses mWindowStatuses;
    private final CobaltCalculation.TargetRegions mEnricher;
    private final CobaltConfig mConfig;

    public CobaltCalculator(final ListMultimap<Chromosome, ReadDepth> mTumorReadDepths, CobaltConfig config) throws IOException
    {
        this.mTumorReadDepths = mTumorReadDepths;
        mWindowStatuses = new WindowStatuses(config);
        mConfig = config;
        mEnricher = mConfig.targetRegionEnricher();
    }

    public ListMultimap<Chromosome, CobaltRatio> doCalculation()
    {
        CobaltCalculation tumorCalculation = new CobaltCalculation(mWindowStatuses, mConfig.RefGenVersion,mEnricher);
        mTumorReadDepths.forEach((tumorCalculation::addReading));
        ListMultimap<Chromosome, CobaltRatio> tumorRatios = tumorCalculation.calculateRatios(mEnricher);
        return tumorRatios;
    }
}
