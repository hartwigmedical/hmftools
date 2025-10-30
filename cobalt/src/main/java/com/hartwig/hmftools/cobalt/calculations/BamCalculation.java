package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_BUCKET_MAX;
import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_BUCKET_MIN;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.cobalt.targeted.CobaltScope;
import com.hartwig.hmftools.common.cobalt.GcMedianReadDepth;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public abstract class BamCalculation
{
    private final ListMultimap<Chromosome, CobaltWindow> WindowsByChromosome = ArrayListMultimap.create();
    private final GCPailsList mGCPailsList = new GCPailsList();
    private final GenomeFilter mGenomeFilter;
    protected final CobaltScope Scope;
    private GcBucketStatistics BucketStatistics;
    final ReadDepthStatisticsNormaliser MeanNormaliser;
    final ResultsNormaliser MegaBaseScaleNormaliser;
    private final ResultsNormaliser FinalNormaliser;

    public BamCalculation(final GenomeFilter mGenomeFilter, CobaltScope scope, RefGenomeVersion version)
    {
        this.mGenomeFilter = mGenomeFilter;
        this.Scope = scope;
        MeanNormaliser = createReadDepthsNormaliser();
        FinalNormaliser = createFinalNormaliser();
        MegaBaseScaleNormaliser = createMegaBaseScaleNormaliser(version);
    }

    public void addReading(Chromosome chromosome, DepthReading readDepth)
    {
        // The genome filter takes into account gc mappability, excluded pseudo-gene regions
        // and excluded non-diploid regions, depending on the mode.
        final boolean isExcluded = mGenomeFilter.exclude(chromosome, readDepth);
        // All windows will be on-target in whole-genome mode.
        final boolean isInTargetRegion = Scope.onTarget(chromosome, readDepth.StartPosition);
        final CobaltWindow rawWindow = new CobaltWindow(chromosome, readDepth, isExcluded, isInTargetRegion);
        final GCPail bucket = mGCPailsList.getGCPail(readDepth.ReadGcContent);
        // Assign a bucket to the window, recording the GC statistics at the same time.
        // The bucket will be set to null in the bucketed window for excluded or off-target readings.
        final CobaltWindow bucketedWindow = rawWindow.bucketed(bucket);
        WindowsByChromosome.put(chromosome, bucketedWindow);
    }

    public ListMultimap<Chromosome, BamRatio> calculateRatios()
    {
        BucketStatistics = new GcBucketStatistics(mGCPailsList, GC_BUCKET_MIN, GC_BUCKET_MAX);
        final ListMultimap<Chromosome, BamRatio> bamResults = ArrayListMultimap.create();
        WindowsByChromosome.forEach((chromosome, window) ->
        {
            // Windows are converted to BamRatios. The essential values of a BamRatio are its ratio
            // and whether it is included (not masked out due to mapability, GC, targeted mode).
            // Subsequent normalisation steps adjust the ratio value, but non-included BamRatios have a ratio fixed at -1.
            BamRatio bamRatio = new BamRatio(chromosome, window.mDepthReading, window.include());
            // Normalise by the GC factor for the window's bucket. This will set the ratio to -1.0 it the bucket is null.
            bamRatio.normaliseForGc(BucketStatistics.medianReadDepth(window.GcBucket));
            // Apply enrichment (does nothing in whole genome mode).
            bamRatio.applyEnrichment(Scope.enrichmentQuotient(chromosome, window.mDepthReading));
            bamResults.put(chromosome, bamRatio);
        });
        // Apply the remaining normalisation and consolidation steps. Some of these are no-ops depending on the mode.
        BamRatios bamRatios = new BamRatios(bamResults);
        bamRatios.normalise(MeanNormaliser);
        bamRatios.consolidate(consolidator());
        bamRatios.normalise(MegaBaseScaleNormaliser);
        bamRatios.normalise(FinalNormaliser);
        return bamRatios.Ratios;
    }

    GcMedianReadDepth medianReadDepths()
    {
        return new GcMedianReadDepth(MeanNormaliser.readDepthMean(), MeanNormaliser.readDepthMedian(), BucketStatistics.bucketToMedianReadDepth());
    }

    abstract ReadDepthStatisticsNormaliser createReadDepthsNormaliser();

    abstract ResultsNormaliser createMegaBaseScaleNormaliser(RefGenomeVersion version);

    ResultsNormaliser createFinalNormaliser()
    {
        return Scope.finalNormaliser();
    }

    abstract ResultsConsolidator consolidator();
}
