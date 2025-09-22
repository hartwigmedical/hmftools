package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_BUCKET_MAX;
import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_BUCKET_MIN;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.cobalt.targeted.TargetRegions;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class CobaltCalculation
{

    private final ListMultimap<Chromosome, CobaltWindow> mWindowsByChromosome = ArrayListMultimap.create();
    private final GCPailsList mGCPailsList = new GCPailsList();
    private final GenomeFilter mGenomeFilter;
    private final RefGenomeVersion GenomeVersion;
    private final TargetRegions mTargetRegions;

    public CobaltCalculation(final GenomeFilter mGenomeFilter, RefGenomeVersion genomeVersion, TargetRegions targetRegions)
    {
        this.mGenomeFilter = mGenomeFilter;
        this.mTargetRegions = targetRegions;
        this.GenomeVersion = genomeVersion;
    }

    public void addReading(Chromosome chromosome, DepthReading readDepth)
    {
        final boolean isExcluded = mGenomeFilter.exclude(chromosome, readDepth);
        final boolean isInTargetRegion = mTargetRegions.onTarget(chromosome, readDepth.StartPosition);
        final CobaltWindow rawWindow = new CobaltWindow(chromosome, readDepth, isExcluded, isInTargetRegion);
        GCPail bucket = mGCPailsList.getGCPail(readDepth.ReadGcContent);
        mWindowsByChromosome.put(chromosome, rawWindow.bucketed(bucket));
    }

    public ListMultimap<Chromosome, CobaltRatio> calculateRatios()
    {
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(mGCPailsList, GC_BUCKET_MIN, GC_BUCKET_MAX);
        ResultsNormaliser finalNormaliser = mTargetRegions.createNormaliser();
        final ListMultimap<Chromosome, BamRatio> bamResults = ArrayListMultimap.create();
        mWindowsByChromosome.forEach((chromosome, window) ->
        {
            BamRatio bamRatio = new BamRatio(chromosome, window.mDepthReading, mTargetRegions.onTarget(chromosome, window.Position));
            bamRatio.normaliseForGc(bucketStatistics.medianReadDepth(window.GcBucket));
            bamRatio.applyEnrichment(mTargetRegions.enrichmentQuotient(chromosome, window.mDepthReading));
            finalNormaliser.recordValue(bamRatio);
            bamResults.put(chromosome, bamRatio);
        });

        final ListMultimap<Chromosome, CobaltRatio> finalResults = ArrayListMultimap.create();
        bamResults.forEach(((chromosome, bamRatio) ->
        {
            finalNormaliser.applyNormalisation(bamRatio);
            finalResults.put(chromosome, bamRatio.toTumorRatio(GenomeVersion));
        }));
        return finalResults;
    }
}
