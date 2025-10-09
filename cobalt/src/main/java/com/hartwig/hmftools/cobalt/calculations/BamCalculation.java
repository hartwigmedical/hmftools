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
        // Use the supplied filters to set the state of the reading and assign it to a GC bucket.
        final boolean isExcluded = mGenomeFilter.exclude(chromosome, readDepth);
        final boolean isInTargetRegion = Scope.onTarget(chromosome, readDepth.StartPosition);
        final GCPail bucket = mGCPailsList.getGCPail(readDepth.ReadGcContent);
        final CobaltWindow rawWindow = new CobaltWindow(chromosome, readDepth, isExcluded, isInTargetRegion);
        WindowsByChromosome.put(chromosome, rawWindow.bucketed(bucket));
    }

    public ListMultimap<Chromosome, BamRatio> calculateRatios()
    {
        BucketStatistics = new GcBucketStatistics(mGCPailsList, GC_BUCKET_MIN, GC_BUCKET_MAX);
        final ListMultimap<Chromosome, BamRatio> bamResults = ArrayListMultimap.create();
        WindowsByChromosome.forEach((chromosome, window) ->
        {
            // Create a reading for each window, normalise according to its GC bucket,
            // then enrich according to supplied per-window values.
            BamRatio bamRatio = new BamRatio(chromosome, window.mDepthReading, Scope.onTarget(chromosome, window.Position));
            bamRatio.normaliseForGc(BucketStatistics.medianReadDepth(window.GcBucket));
            bamRatio.applyEnrichment(Scope.enrichmentQuotient(chromosome, window.mDepthReading));
            bamResults.put(chromosome, bamRatio);
        });
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
