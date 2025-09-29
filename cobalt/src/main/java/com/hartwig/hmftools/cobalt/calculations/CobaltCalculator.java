package com.hartwig.hmftools.cobalt.calculations;

import java.util.List;

import com.google.common.base.Preconditions;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.CobaltConfig;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.cobalt.targeted.CobaltScope;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class CobaltCalculator
{
    private final ListMultimap<Chromosome, CobaltRatio> mRatios;
    private final List<MedianRatio> mMedianRatios;
    public CobaltCalculator(
            final ListMultimap<Chromosome, DepthReading> tumourDepthReadings,
            final ListMultimap<Chromosome, DepthReading> referenceDepthReadings,
            CobaltConfig config)
    {
        Preconditions.checkArgument(!tumourDepthReadings.isEmpty() || !referenceDepthReadings.isEmpty());
        GenomeFilter mWindowStatuses = new WindowStatuses(config.gcProfileData(), config.excludedRegions());
        CobaltScope enricher = config.scope();

        TumorBamCalculation tumorBamCalculation = new TumorBamCalculation(mWindowStatuses, enricher);
        tumourDepthReadings.forEach((tumorBamCalculation::addReading));
        ListMultimap<Chromosome, BamRatio> tumorResults = tumorBamCalculation.calculateRatios();

        ReferenceBamCalculation referenceBamCalculation = new ReferenceBamCalculation(mWindowStatuses, enricher);
        referenceDepthReadings.forEach((referenceBamCalculation::addReading));
        ListMultimap<Chromosome, BamRatio> referenceResults = referenceBamCalculation.calculateRatios();
        mMedianRatios = referenceBamCalculation.medianRatios(config.genomeVersion());

        ResultsCollator collator = new ResultsCollator(config.genomeVersion());
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
}
