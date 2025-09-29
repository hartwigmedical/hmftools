package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_BUCKET_MAX;
import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_BUCKET_MIN;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.cobalt.targeted.CobaltScope;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public abstract class BamCalculation
{
    private final ListMultimap<Chromosome, CobaltWindow> mWindowsByChromosome = ArrayListMultimap.create();
    private final GCPailsList mGCPailsList = new GCPailsList();
    private final GenomeFilter mGenomeFilter;
    protected final CobaltScope mScope;

    public BamCalculation(final GenomeFilter mGenomeFilter, CobaltScope scope)
    {
        this.mGenomeFilter = mGenomeFilter;
        this.mScope = scope;
    }

    public void addReading(Chromosome chromosome, DepthReading readDepth)
    {
        final boolean isExcluded = mGenomeFilter.exclude(chromosome, readDepth);
        final boolean isInTargetRegion = mScope.onTarget(chromosome, readDepth.StartPosition);
        final GCPail bucket = mGCPailsList.getGCPail(readDepth.ReadGcContent);
        final CobaltWindow rawWindow = new CobaltWindow(chromosome, readDepth, isExcluded, isInTargetRegion);
        mWindowsByChromosome.put(chromosome, rawWindow.bucketed(bucket));
    }

    public ListMultimap<Chromosome, BamRatio> calculateRatios()
    {
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(mGCPailsList, GC_BUCKET_MIN, GC_BUCKET_MAX);
        ResultsNormaliser finalNormaliser = finalMeanNormaliser();
        ResultsNormaliser diploidNormaliser = diploidNormaliser();
        final ListMultimap<Chromosome, BamRatio> bamResults = ArrayListMultimap.create();
        mWindowsByChromosome.forEach((chromosome, window) ->
        {
            BamRatio bamRatio = new BamRatio(chromosome, window.mDepthReading, mScope.onTarget(chromosome, window.Position));
            bamRatio.normaliseForGc(bucketStatistics.medianReadDepth(window.GcBucket));
            bamRatio.applyEnrichment(mScope.enrichmentQuotient(chromosome, window.mDepthReading));
            diploidNormaliser.recordValue(bamRatio);
            finalNormaliser.recordValue(bamRatio);
            bamResults.put(chromosome, bamRatio);
        });

        diploidNormaliser.recordsAllAdded();
        finalNormaliser.recordsAllAdded();

        bamResults.forEach(((chromosome, bamRatio) ->
        {
            diploidNormaliser.applyNormalisation(bamRatio);
            finalNormaliser.applyNormalisation(bamRatio);
        }));
        return bamResults;
    }

    abstract ResultsNormaliser finalMeanNormaliser();

    abstract ResultsNormaliser diploidNormaliser();
}
