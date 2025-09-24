package com.hartwig.hmftools.cobalt.calculations;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.CobaltConfig;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.cobalt.targeted.TargetRegions;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class CobaltCalculator
{
    private final ListMultimap<Chromosome, CobaltRatio> mRatios;

    public CobaltCalculator(
            final ListMultimap<Chromosome, DepthReading> tumourDepthReadings,
            final ListMultimap<Chromosome, DepthReading> referenceDepthReadings,
            CobaltConfig config)
    {
        GenomeFilter mWindowStatuses = new WindowStatuses(config.gcProfileData(), config.mExcludedRegions);
        TargetRegions mEnricher = config.targetRegionEnricher();

        TumorBamCalculation tumorBamCalculation = new TumorBamCalculation(mWindowStatuses, mEnricher);
        tumourDepthReadings.forEach((tumorBamCalculation::addReading));
        ListMultimap<Chromosome, BamRatio> tumorResults = tumorBamCalculation.calculateRatios();

        ReferenceBamCalculation referenceBamCalculation = new ReferenceBamCalculation(mWindowStatuses, mEnricher);
        referenceDepthReadings.forEach((referenceBamCalculation::addReading));
        ListMultimap<Chromosome, BamRatio> referenceResults = referenceBamCalculation.calculateRatios();

        ResultsCollator collator = new ResultsCollator(config.RefGenVersion);
        mRatios = collator.collateResults(tumorResults, referenceResults);
    }

    public ListMultimap<Chromosome, CobaltRatio> getCalculatedRatios()
    {
        return mRatios;
    }
}
