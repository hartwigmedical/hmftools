package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_BUCKET_MAX;
import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_BUCKET_MIN;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.CobaltUtils;
import com.hartwig.hmftools.cobalt.count.ReadDepth;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

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

        boolean applyFinalNormalisation();
    }

    private final ListMultimap<Chromosome, CobaltWindow> mWindowsByChromosome = ArrayListMultimap.create();
    private final GCPailsList mGCPailsList = new GCPailsList();
    private final Filter mFilter;

    public CobaltCalculation(final Filter mFilter)
    {
        this.mFilter = mFilter;
    }

    public void addReading(Chromosome chromosome, ReadDepth readDepth)
    {
        //        Preconditions.checkArgument(chromosome.equals(readDepth.Chromosome)); // todo reinstate after making ReadDepth have a Chromosome
        if(mFilter.exclude(chromosome, readDepth))
        {
            return;
        }
        GCPail bucket = mGCPailsList.getGCPail(readDepth.ReadGcContent);
        bucket.addReading(readDepth.ReadDepth);
        CobaltWindow window = new CobaltWindow(chromosome, readDepth.StartPosition, readDepth, bucket);
        mWindowsByChromosome.put(chromosome, window);
    }

    public ListMultimap<Chromosome, CobaltRatio> calculateRatios(TargetRegions targetRegions)
    {
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(mGCPailsList, GC_BUCKET_MIN, GC_BUCKET_MAX);
        DescriptiveStatistics meanCalculator = new DescriptiveStatistics();
        final ListMultimap<Chromosome, CobaltRatio> interimResults = ArrayListMultimap.create();
        mWindowsByChromosome.forEach((chromosome, window) ->
        {
            // Maybe use something other tnan CobaltRatio here - something that can
            // be transformed in a logical way by the arithmetical normalisation operations
            // including handling getting set to NaN or whatever
            if(targetRegions.isInTargetRegions(chromosome, window))
            {
                double medianReadDepthForBucket = bucketStatistics.medianReadDepth(window.GcBucket.mGC);
                if (medianReadDepthForBucket < 0)
                {
                    // disallowed bucket
                    // todo - better just make normalise(-1.0)  do this
                    interimResults.put(chromosome, CobaltUtils.tumorOnlyRatio(window, -1.0));
                }
                else
                {
                    double ratio = window.ReadDepth.ReadDepth / medianReadDepthForBucket;
                    ratio = ratio / targetRegions.enrichmentQuotient(chromosome, window.ReadDepth);
                    if(targetRegions.applyFinalNormalisation())
                    {
                        meanCalculator.addValue(ratio);
                    }
                    interimResults.put(chromosome, CobaltUtils.tumorOnlyRatio(window, ratio));
                }
            }
        });
        ListMultimap<Chromosome, CobaltRatio> finalResults;
        if(targetRegions.applyFinalNormalisation())
        {
            final ListMultimap<Chromosome, CobaltRatio> normalisedResults = ArrayListMultimap.create();
            double meanRatio = meanCalculator.getMean();
            interimResults.forEach(((chromosome, cobaltRatio) ->
            {
                if (cobaltRatio.tumorGCRatio() < 0)
                {
                    normalisedResults.put(chromosome, cobaltRatio); // todo - better just make normalise(-1.0)  do nothing
                }
                else
                {
                    CobaltRatio normalised = cobaltRatio.normaliseTumorGcRatio(meanRatio);
                    normalisedResults.put(chromosome, normalised);
                }
            }));
            finalResults = normalisedResults;
        }
        else
        {
            finalResults = interimResults;
        }
        return finalResults;
    }
}
