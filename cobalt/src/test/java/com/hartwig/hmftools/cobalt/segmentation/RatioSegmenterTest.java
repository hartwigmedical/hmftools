package com.hartwig.hmftools.cobalt.segmentation;

import static com.hartwig.hmftools.cobalt.segmentation.Arm.P;
import static com.hartwig.hmftools.cobalt.segmentation.Arm.Q;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.segmentation.PiecewiseConstantFit;

import org.junit.Assert;
import org.junit.Test;

public class RatioSegmenterTest
{
    @Test
    public void segmentationTest() throws Exception
    {
        List<CobaltRatio> ratios = new ArrayList<>();
        // Since we're not testing the segmenter itself,
        // just use the same ratios for each chromosome.
        List<HumanChromosome> chromosomes = List.of(_1, _2, _3);
        chromosomes.forEach(chromosome -> {
            ratios.add(cr(chromosome, 1, 0.5));
            ratios.add(cr(chromosome, 1001, 0.51));
            ratios.add(cr(chromosome, 2001, 0.49));
            ratios.add(cr(chromosome, 3001, -1.0));
            ratios.add(cr(chromosome, 4001, -1.0));
            ratios.add(cr(chromosome, 5001, -1.0));
            ratios.add(cr(chromosome, 6001, -1.0));
            ratios.add(cr(chromosome, 13001, 10.09));
            ratios.add(cr(chromosome, 14001, 9.92));
            ratios.add(cr(chromosome, 15001, 10.02));
            ratios.add(cr(chromosome, 16001, -1.0));

            ratios.add(cr(chromosome, 20_001, 0.5));
            ratios.add(cr(chromosome, 21_001, 0.51));
            ratios.add(cr(chromosome, 22_001, 0.49));
            ratios.add(cr(chromosome, 23_001, 0.52));
        });

        ChrArmLocator locator = cobaltRatio ->
        {
            Arm arm = cobaltRatio.position() > 20_000 ? Q : P;
            return new ChrArm(cobaltRatio.chr(), arm);
        };
        ExecutorService executor = Executors.newFixedThreadPool(2);
        RatioSegmenter segmenter = new RatioSegmenter(ratios, locator, 50.0);
        Map<ChrArm, PiecewiseConstantFit> segments = segmenter.getSegmentation(executor);
        Assert.assertEquals(6, segments.size());
        final PiecewiseConstantFit pFit = new PiecewiseConstantFit(new int[]{3,3}, new int[]{0,3}, new double[]{0.5, 10.01} );
        Assert.assertEquals(pFit, segments.get(ca(_1, P)));
        Assert.assertEquals(pFit, segments.get(ca(_2, P)));
        Assert.assertEquals(pFit, segments.get(ca(_3, P)));
        final PiecewiseConstantFit qFit = new PiecewiseConstantFit(new int[]{4}, new int[]{0}, new double[]{0.505} );
        Assert.assertEquals(qFit, segments.get(ca(_1, Q)));
        Assert.assertEquals(qFit, segments.get(ca(_2, Q)));
        Assert.assertEquals(qFit, segments.get(ca(_3, Q)));
    }

    private static ChrArm ca(HumanChromosome chromosome, Arm arm)
    {
        return new ChrArm(chromosome, arm);
    }
    private static CobaltRatio cr(HumanChromosome chromosome, int start, double value)
    {
        // The RatioSegmenter uses the log of the values.
        // So that we know how the segmentation works out, we will anti-log them first.
        double v = value >= 0 ? Math.pow(2, value) : -1.0;
        return new CobaltRatio(chromosome.shortName(), start, v, v, v, v, v, v, v);
    }
}
