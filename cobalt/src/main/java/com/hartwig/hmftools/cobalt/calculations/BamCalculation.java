package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_BUCKET_MAX;
import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_BUCKET_MIN;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.cobalt.targeted.CobaltScope;
import com.hartwig.hmftools.common.cobalt.GcMedianReadDepth;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public abstract class BamCalculation
{
    private final ListMultimap<Chromosome, CobaltWindow> mWindowsByChromosome = ArrayListMultimap.create();
    private final GCPailsList mGCPailsList = new GCPailsList();
    private final GenomeFilter mGenomeFilter;
    protected final CobaltScope mScope;
    private final GcBucketStatistics bucketStatistics = new GcBucketStatistics(mGCPailsList, GC_BUCKET_MIN, GC_BUCKET_MAX);
    private final ReadDepthStatisticsNormaliser meanNormaliser = createReadDepthsNormaliser();
    private final ResultsNormaliser megaBaseScaleNormaliser = createMegaBaseScaleNormaliser();
    private final ResultsNormaliser finalNormaliser = createFinalNormaliser();

    public BamCalculation(final GenomeFilter mGenomeFilter, CobaltScope scope)
    {
        this.mGenomeFilter = mGenomeFilter;
        this.mScope = scope;
    }

    public void addReading(Chromosome chromosome, DepthReading readDepth)
    {
        // Use the supplied filters to set the state of the reading and assign it to a GC bucket.
        final boolean isExcluded = mGenomeFilter.exclude(chromosome, readDepth);
        final boolean isInTargetRegion = mScope.onTarget(chromosome, readDepth.StartPosition);
        final GCPail bucket = mGCPailsList.getGCPail(readDepth.ReadGcContent);
        final CobaltWindow rawWindow = new CobaltWindow(chromosome, readDepth, isExcluded, isInTargetRegion);
        mWindowsByChromosome.put(chromosome, rawWindow.bucketed(bucket));
    }

    public ListMultimap<Chromosome, BamRatio> calculateRatios()
    {
        final ListMultimap<Chromosome, BamRatio> bamResults = ArrayListMultimap.create();
        mWindowsByChromosome.forEach((chromosome, window) ->
        {
            // Normalise each reading according to its GC bucket, then enrich according
            // to supplied per-window values, then record this result in the mean normaliser.
            BamRatio bamRatio = new BamRatio(chromosome, window.mDepthReading, mScope.onTarget(chromosome, window.Position));
            bamRatio.normaliseForGc(bucketStatistics.medianReadDepth(window.GcBucket));
            bamRatio.applyEnrichment(mScope.enrichmentQuotient(chromosome, window.mDepthReading));
            meanNormaliser.recordValue(bamRatio);
            bamResults.put(chromosome, bamRatio);
        });

        // Apply the mean normaliser and record the normalised values in the mega-base scale normaliser.
        meanNormaliser.dataCollectionFinished();
        bamResults.forEach(((chromosome, bamRatio) ->
        {
            meanNormaliser.normalise(bamRatio);
            megaBaseScaleNormaliser.recordValue(bamRatio);
        }));

        // Apply the mega-base scale normaliser and record the normalised value for final normalisation.
        megaBaseScaleNormaliser.dataCollectionFinished();
        bamResults.forEach(((chromosome, bamRatio) ->
        {
            megaBaseScaleNormaliser.normalise(bamRatio);
            finalNormaliser.recordValue(bamRatio);
        }));

        // Apply the final normalisation step.
        finalNormaliser.dataCollectionFinished();
        bamResults.forEach(((chromosome, bamRatio) -> finalNormaliser.normalise(bamRatio)));
        return bamResults;
    }

    GcMedianReadDepth medianReadDepths()
    {
        // todo test
        return new GcMedianReadDepth(meanNormaliser.readDepthMean(), meanNormaliser.readDepthMedian(), mGCPailsList.bucketToMedianReadDepth());
    }

    abstract ReadDepthStatisticsNormaliser createReadDepthsNormaliser();

    abstract ResultsNormaliser createMegaBaseScaleNormaliser();

    abstract ResultsNormaliser createFinalNormaliser();
}
