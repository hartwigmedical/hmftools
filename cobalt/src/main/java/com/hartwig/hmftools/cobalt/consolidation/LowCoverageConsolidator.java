package com.hartwig.hmftools.cobalt.consolidation;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator.consolidateIntoBuckets;

import java.util.Iterator;
import java.util.List;

import com.google.common.base.Preconditions;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.calculations.BamRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class LowCoverageConsolidator implements ResultsConsolidator
{
    private final ArrayListMultimap<Chromosome, LowCovBucket> Boundaries = ArrayListMultimap.create();
    private final int ConsolidationCount;

    public LowCoverageConsolidator(final int consolidationCount)
    {
        Preconditions.checkArgument(consolidationCount > 1);
        ConsolidationCount = consolidationCount;
    }

    @Override
    public ListMultimap<Chromosome, BamRatio> consolidate(ListMultimap<Chromosome, BamRatio> ratios)
    {
        createBoundariesIfNecessary(ratios);
        return populateLowCoverageRatio(ratios);
    }

    private void createBoundariesIfNecessary(ListMultimap<Chromosome, BamRatio> ratios)
    {
        if(!Boundaries.isEmpty())
        {
            return;
        }
        for(Chromosome chromosome : ratios.keySet())
        {
            List<Integer> nonMaskedPositions = ratios.get(chromosome).stream()
                    .filter(bamRatio -> bamRatio.ratio() >= 0)
                    .map(bamRatio -> bamRatio.Position).toList();
            List<LowCovBucket> consolidatedBuckets = consolidateIntoBuckets(nonMaskedPositions, ConsolidationCount);
            Boundaries.putAll(chromosome, consolidatedBuckets);
            CB_LOGGER.info("chromosome: {}, low cov buckets count: {}", chromosome, consolidatedBuckets.size());
        }
    }

    // we create a pan window ratio by taking the mean count of super windows that combine multiple windows
    private ListMultimap<Chromosome, BamRatio> populateLowCoverageRatio(final ListMultimap<Chromosome, BamRatio> rawRatios)
    {
        ListMultimap<Chromosome, BamRatio> result = ArrayListMultimap.create();
        for(Chromosome chromosome : Boundaries.keySet())
        {
            List<BamRatio> chromosomeRatios = rawRatios.get(chromosome);
            Iterator<LowCovBucket> bucketItr = Boundaries.get(chromosome).iterator();
            LowCovBucket bucket;
            DescriptiveStatistics bucketRatios = new DescriptiveStatistics();
            DescriptiveStatistics bucketGCs = new DescriptiveStatistics();
            if(!bucketItr.hasNext())
            {
                CB_LOGGER.error("low cov bucket for chromosome {} not found", chromosome);
                continue;
            }
            else
            {
                bucket = bucketItr.next();
            }

            for(BamRatio bamRatio : chromosomeRatios)
            {
                if(bucket == null)
                {
                    // no bucket for whole chromosome, or we already finished last bucket
                    continue;
                }
                if(bamRatio.ratio() >= 0)
                {
                    if(bamRatio.Position > bucket.EndPosition)
                    {
                        // Turn the current bucket into a consolidated ratio and add it to the results
                        BamRatio consolidated =
                                new BamRatio(chromosome, bucket.BucketPosition, bucketRatios.getMean(), bucketGCs.getMean());
                        result.put(chromosome, consolidated);
                        if(bucketItr.hasNext())
                        {
                            // move to next bucket
                            bucket = bucketItr.next();
                            bucketRatios = new DescriptiveStatistics();
                            bucketGCs = new DescriptiveStatistics();
                            bucketRatios.addValue(bamRatio.ratio());
                            bucketGCs.addValue(bamRatio.gcContent());
                        }
                        else
                        {
                            // no more bucket for this chromosome. Setting bucket to null will let us skip through
                            // the rest of the chromosome
                            bucket = null;
                        }
                    }
                    else
                    {
                        // Add the data for this ratio to the current bucket data.
                        bucketRatios.addValue(bamRatio.ratio());
                        bucketGCs.addValue(bamRatio.gcContent());
                    }
                }
            }
            // Wrap up the last bucket, if there is one.
            if(bucket != null)
            {
                BamRatio consolidated = new BamRatio(chromosome, bucket.BucketPosition, bucketRatios.getMean(), bucketGCs.getMean());
                result.put(chromosome, consolidated);
            }
        }
        return result;
    }
}
