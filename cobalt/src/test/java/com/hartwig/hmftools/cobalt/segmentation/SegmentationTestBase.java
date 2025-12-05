package com.hartwig.hmftools.cobalt.segmentation;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._X;
import static com.hartwig.hmftools.common.segmentation.Arm.P;
import static com.hartwig.hmftools.common.segmentation.Arm.Q;

import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.segmentation.Arm;
import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

import org.junit.After;
import org.junit.Before;

public class SegmentationTestBase
{
    ListMultimap<Chromosome, CobaltRatio> ratios;
    ExecutorService executor;
    ChrArmLocator locator = cobaltRatio ->
    {
        Arm arm = cobaltRatio.position() > 20_000 ? Q : P;
        return new ChrArm(cobaltRatio.chr(), arm);
    };

    @Before
    public void setup()
    {
        executor = Executors.newFixedThreadPool(8);
        List<HumanChromosome> chromosomes = List.of(_1, _2, _X);
        ratios = ArrayListMultimap.create();
        chromosomes.forEach(chromosome ->
        {
            ratios.put(chromosome, ratioForValue(chromosome, 1, 0.5));
            ratios.put(chromosome, ratioForValue(chromosome, 1001, 0.51));
            ratios.put(chromosome, ratioForValue(chromosome, 2001, 0.49));
        });
    }

    @After
    public void cleanup()
    {
        executor.shutdown();
    }

    CobaltRatio ratioForValue(HumanChromosome chromosome, int start, double value)
    {
        // The RatioSegmenter uses the log of the values.
        // So that we know how the segmentation works out, we will anti-log them first.
        double v = value >= 0 ? Math.pow(2, value) : -1.0;
        return new CobaltRatio(chromosome.shortName(), start, v, v, v, v, v, v, v);
    }

    CobaltRatio ratio(HumanChromosome chromosome, int start, double v)
    {
        return new CobaltRatio(chromosome.shortName(), start, v, v, v, v, v, v, v);
    }
}
