package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_BUCKET_MAX;
import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_BUCKET_MIN;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.count.ReadDepth;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class CobaltCalculation
{
    public interface Filter
    {
        boolean exclude(final Chromosome chromosome, ReadDepth readDepth);
    }

    public interface TargetRegions
    {
        boolean isInTargetRegions(final Chromosome chromosome, CobaltWindow window);

        double enrichmentQuotient(Chromosome chromosome, ReadDepth readDepth);

        boolean isInTargetRegions(Chromosome chromosome, int position);

        boolean applyFinalNormalisation();
    }

    private final ListMultimap<Chromosome, CobaltWindow> mWindowsByChromosome = ArrayListMultimap.create();
    private final GCPailsList mGCPailsList = new GCPailsList();
    private final Filter mFilter;
    private final RefGenomeVersion GenomeVersion;
    private final TargetRegions mTargetRegions;

    public CobaltCalculation(final Filter mFilter, RefGenomeVersion genomeVersion, TargetRegions targetRegions)
    {
        this.mFilter = mFilter;
        this.mTargetRegions = targetRegions;
        this.GenomeVersion = genomeVersion;
    }

    public void addReading(Chromosome chromosome, ReadDepth readDepth)
    {
        CobaltWindow window;

        if(mFilter.exclude(chromosome, readDepth))
        {
            window = new CobaltWindow(chromosome, readDepth.StartPosition, readDepth);
        }
        else
        {
            GCPail bucket = mGCPailsList.getGCPail(readDepth.ReadGcContent);
            if(chromosome.isAutosome() && mTargetRegions.isInTargetRegions(chromosome, readDepth.StartPosition) && readDepth.ReadDepth > 0)
            {
                bucket.addReading(readDepth.ReadDepth);
            }
            window = new CobaltWindow(chromosome, readDepth.StartPosition, readDepth, bucket);
        }
        mWindowsByChromosome.put(chromosome, window);
    }

    public ListMultimap<Chromosome, CobaltRatio> calculateRatios(TargetRegions targetRegions)
    {
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(mGCPailsList, GC_BUCKET_MIN, GC_BUCKET_MAX);
        BamRatiosMeanCalculator meanCalculator = new BamRatiosMeanCalculator();
        final ListMultimap<Chromosome, BamRatio> bamResults = ArrayListMultimap.create();
        mWindowsByChromosome.forEach((chromosome, window) ->
        {
            BamRatio bamRatio = new BamRatio(chromosome, window.ReadDepth, targetRegions.isInTargetRegions(chromosome, window));
            bamRatio.normaliseForGc(bucketStatistics.medianReadDepth(window.GcBucket));
            bamRatio.applyEnrichment(targetRegions.enrichmentQuotient(chromosome, window.ReadDepth));

            if(targetRegions.applyFinalNormalisation())
            {
                meanCalculator.recordValue(bamRatio);
            }
            bamResults.put(chromosome, bamRatio);
        });

        final ListMultimap<Chromosome, CobaltRatio> finalResults = ArrayListMultimap.create();
        if(targetRegions.applyFinalNormalisation())
        {
            double meanRatio = meanCalculator.mean();
            CB_LOGGER.info("Normalising results with mean ratio {}", meanRatio);

            bamResults.forEach(((chromosome, cobaltRatio) ->
            {
                cobaltRatio.normaliseByMean(meanRatio);
                finalResults.put(chromosome, cobaltRatio.toTumorRatio(GenomeVersion));
            }));
        }
        else
        {
            bamResults.forEach(((chromosome, cobaltRatio) ->
            {
                finalResults.put(chromosome, cobaltRatio.toTumorRatio(GenomeVersion));
            }));
        }
        return finalResults;
    }
}
