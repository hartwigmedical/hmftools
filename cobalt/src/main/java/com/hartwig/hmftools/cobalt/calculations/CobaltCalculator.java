package com.hartwig.hmftools.cobalt.calculations;

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
    private final GcMedianReadDepth mTumorMedianReadDepths;
    private final GcMedianReadDepth mReferenceMedianReadDepths;
    public CobaltCalculator(
            final ListMultimap<HumanChromosome, DepthReading> tumourDepthReadings,
            final ListMultimap<HumanChromosome, DepthReading> referenceDepthReadings,
            CobaltConfig config)
    {
        Preconditions.checkArgument(!tumourDepthReadings.isEmpty() || !referenceDepthReadings.isEmpty());
        GenomeFilter mWindowStatuses = new WindowStatuses(config.gcProfileData(), config.excludedRegions());
        CobaltScope scope = config.scope();

        TumorBamCalculation tumorBamCalculation = new TumorBamCalculation(mWindowStatuses, scope, config.genomeVersion());
        tumourDepthReadings.forEach(tumorBamCalculation::addReading);
        ListMultimap<Chromosome, BamRatio> tumorResults = tumorBamCalculation.calculateRatios();
        mTumorMedianReadDepths = tumorBamCalculation.medianReadDepths();

        ReferenceBamCalculation referenceBamCalculation = new ReferenceBamCalculation(mWindowStatuses, scope, config.genomeVersion());
        referenceDepthReadings.forEach(referenceBamCalculation::addReading);
        ListMultimap<Chromosome, BamRatio> referenceResults = referenceBamCalculation.calculateRatios();
        mMedianRatios = referenceBamCalculation.medianRatios();
        mReferenceMedianReadDepths = referenceBamCalculation.medianReadDepths();

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

    public GcMedianReadDepth tumorMedianReadDepth()
    {
        return mTumorMedianReadDepths;
    }

    public GcMedianReadDepth referenceMedianReadDepth()
    {
        return mReferenceMedianReadDepths;
    }
}
