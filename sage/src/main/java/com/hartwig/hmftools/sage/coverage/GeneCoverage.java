package com.hartwig.hmftools.sage.coverage;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.apache.commons.compress.utils.Lists;
import org.apache.commons.math3.distribution.PoissonDistribution;

public class GeneCoverage implements Consumer<GenomeRegion>
{
    private final String mChromosome;
    private final String mGene;
    private final List<ExonCoverage> mExonCoverage;
    private final long mMinPosition;
    private final long mMaxPosition;

    public static final List<Integer> DEPTH_BUCKETS = Lists.newArrayList();
    public static final int MAX_DEPTH_BUCKET = 10000;

        /* frequency distribution of depth in units of:
        - 1 up to 30
        - 5 up to 100
        - 50 up to 500  (8 regions)
        - 100 up to 2000 (15 regions)
        - 1000 up to 10000 (8 regions)
        - 10000+ (1 region)
    */

    public static void populateCoverageBuckets()
    {
        if(!DEPTH_BUCKETS.isEmpty())
            return;

        int depthBucket = 0;
        int increment = 1;
        for(; depthBucket < 30; depthBucket += increment)
        {
            DEPTH_BUCKETS.add(depthBucket);
        }

        increment = 10;
        for(; depthBucket < 100; depthBucket += increment)
        {
            DEPTH_BUCKETS.add(depthBucket);
        }

        increment = 50;
        for(; depthBucket < 500; depthBucket += increment)
        {
            DEPTH_BUCKETS.add(depthBucket);
        }

        increment = 100;
        for(; depthBucket < 2000; depthBucket += increment)
        {
            DEPTH_BUCKETS.add(depthBucket);
        }

        increment = 1000;
        for(; depthBucket <= 10000; depthBucket += increment)
        {
            DEPTH_BUCKETS.add(depthBucket);
        }
    }

    public GeneCoverage(final List<NamedBed> exons)
    {
        mChromosome = exons.get(0).chromosome();
        mGene = exons.get(0).name();
        mExonCoverage = exons.stream().map(ExonCoverage::new).collect(Collectors.toList());

        long tmpMin = exons.get(0).start();
        long tmpMax = exons.get(0).end();

        for(NamedBed exon : exons)
        {
            tmpMin = Math.min(tmpMin, exon.start());
            tmpMax = Math.max(tmpMax, exon.end());
        }

        mMinPosition = tmpMin;
        mMaxPosition = tmpMax;
    }

    @Override
    public void accept(final GenomeRegion alignment)
    {
        if(alignment.chromosome().equals(mChromosome))
        {
            if(alignment.start() <= mMaxPosition && alignment.end() >= mMinPosition)
            {
                mExonCoverage.forEach(x -> x.accept(alignment));
            }
        }
    }

    public String chromosome()
    {
        return mChromosome;
    }

    public GeneDepth geneDepth()
    {
        int[] depthCounts = baseCoverageSummary(mExonCoverage);

        return new GeneDepth(mGene, missedVariantLikelihood(depthCounts), depthCounts);
    }

    private static int[] baseCoverageSummary(final Collection<ExonCoverage> exons)
    {
        int[] geneDepth = new int[DEPTH_BUCKETS.size()];

        for(ExonCoverage exon : exons)
        {
            for(int baseDepth : exon.coverage())
            {
                geneDepth[bucket(baseDepth)]++;
            }
        }

        return geneDepth;
    }

    public static int depth(int bucket)
    {
        if(bucket >= DEPTH_BUCKETS.size())
            return MAX_DEPTH_BUCKET;

        if(bucket == DEPTH_BUCKETS.size() - 1)
            return DEPTH_BUCKETS.get(DEPTH_BUCKETS.size() - 1);

        int depth = DEPTH_BUCKETS.get(bucket);
        int depthNext = DEPTH_BUCKETS.get(bucket + 1);

        if(depthNext == depth + 1)
            return depth;

        return (depth + depthNext) / 2;
    }

    public static int bucket(int depth)
    {
        for(int i = 0; i < DEPTH_BUCKETS.size(); ++i)
        {
            int depthBucket = DEPTH_BUCKETS.get(i);

            if(depth <= depthBucket)
                return depth < depthBucket ? i - 1 : i;
        }

        return DEPTH_BUCKETS.size() - 1;
    }

    public static double missedVariantLikelihood(int[] baseCoverage)
    {
        int totalCoverage = Arrays.stream(baseCoverage).sum();
        double totalLikelihood = 0;

        for(int i = 0; i < baseCoverage.length; i++)
        {
            int depth = depth(i);
            int coverage = baseCoverage[i];

            if(coverage > 0)
            {
                final double proportion = 1d * coverage / totalCoverage;
                final double likelihoodOfMissing;
                if(depth == 0)
                {
                    likelihoodOfMissing = 1;
                }
                else
                {
                    final PoissonDistribution distribution = new PoissonDistribution(depth / 2d);
                    likelihoodOfMissing = distribution.cumulativeProbability(2);
                }

                totalLikelihood += proportion * likelihoodOfMissing;
            }
        }

        return totalLikelihood;
    }
}
