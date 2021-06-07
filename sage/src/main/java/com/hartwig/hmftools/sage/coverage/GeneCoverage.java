package com.hartwig.hmftools.sage.coverage;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.jetbrains.annotations.NotNull;

public class GeneCoverage implements Consumer<GenomeRegion>
{

    private static final int MAX_BUCKET = 37;

    private final String chromosome;
    private final String gene;
    private final List<ExonCoverage> exonCoverage;
    private final long minPosition;
    private final long maxPosition;

    public GeneCoverage(final List<NamedBed> exons)
    {
        assert (!exons.isEmpty());
        this.chromosome = exons.get(0).chromosome();
        this.gene = exons.get(0).name();
        this.exonCoverage = exons.stream().map(ExonCoverage::new).collect(Collectors.toList());

        long tmpMin = exons.get(0).start();
        long tmpMax = exons.get(0).end();

        for(NamedBed exon : exons)
        {
            tmpMin = Math.min(tmpMin, exon.start());
            tmpMax = Math.max(tmpMax, exon.end());
        }

        minPosition = tmpMin;
        maxPosition = tmpMax;
    }

    @Override
    public void accept(final GenomeRegion alignment)
    {
        if(alignment.chromosome().equals(chromosome))
        {
            if(alignment.start() <= maxPosition && alignment.end() >= minPosition)
            {
                exonCoverage.forEach(x -> x.accept(alignment));
            }
        }
    }

    public String chromosome()
    {
        return chromosome;
    }

    @NotNull
    public GeneDepth geneDepth()
    {
        int[] depthCounts = baseCoverageSummary(exonCoverage);

        return ImmutableGeneDepth.builder()
                .gene(gene)
                .depthCounts(depthCounts)
                .missedVariantLikelihood(missedVariantLikelihood(depthCounts))
                .build();
    }

    static int[] baseCoverageSummary(@NotNull Collection<ExonCoverage> exons)
    {
        int[] geneDepth = new int[MAX_BUCKET + 1];
        for(ExonCoverage exon : exons)
        {
            for(int baseDepth : exon.coverage())
            {
                geneDepth[bucket(baseDepth)]++;
            }
        }

        return geneDepth;
    }

    static double missedVariantLikelihood(int[] baseCoverage)
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

    static int depth(int bucket)
    {
        if(bucket < 30)
        {
            return bucket;
        }

        if(bucket <= 36)
        {
            return (bucket - 30) * 10 + 35;
        }

        return 100;
    }

    static int bucket(int depth)
    {
        if(depth <= 30)
        {
            return depth;
        }

        for(int i = 0; i < 7; i++)
        {
            int maxDepth = 40 + 10 * i;
            if(depth < maxDepth)
            {
                return 30 + i;
            }
        }

        return MAX_BUCKET;
    }
}
