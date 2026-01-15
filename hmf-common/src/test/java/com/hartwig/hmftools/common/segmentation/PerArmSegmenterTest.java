package com.hartwig.hmftools.common.segmentation;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._X;
import static com.hartwig.hmftools.common.segmentation.Arm.P;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class PerArmSegmenterTest
{
    private static final double GAMMA_100 = 100.0;
    ListMultimap<Chromosome, PositionDepth> ratios;
    ExecutorService executor;
    ChrArmLocator locator = cobaltRatio -> new ChrArm(cobaltRatio.chr(), P);

    @Before
    public void setup()
    {
        executor = Executors.newFixedThreadPool(8);
        List<HumanChromosome> chromosomes = List.of(_1, _2, _X);
        ratios = ArrayListMultimap.create();
        chromosomes.forEach(chromosome ->
        {
            ratios.put(chromosome, ratio(chromosome, 1, 0.5));
            ratios.put(chromosome, ratio(chromosome, 1001, 0.51));
            ratios.put(chromosome, ratio(chromosome, 2001, 0.49));
        });
    }

    @Test
    public void uniformPenaltyThresholdTest() throws Exception
    {
        int threshold = 100;
        loadCutoffSensitiveData(threshold - 15);
        PAS segmenter = new PAS(ratios, locator, GAMMA_100, threshold);
        Map<ChrArm, ChromosomeArmSegments<PositionDepth>> armToSegments = segmenter.getSegmentation(executor);
        ChromosomeArmSegments<PositionDepth> segments = armToSegments.get(new ChrArm(_2, P));
        Assert.assertEquals(3, segments.Segments.size());

        loadCutoffSensitiveData(threshold - 14);
        segmenter = new PAS(ratios, locator, GAMMA_100, threshold - 14);
        armToSegments = segmenter.getSegmentation(executor);
        segments = armToSegments.get(new ChrArm(_2, P));
        Assert.assertEquals(4, segments.Segments.size());
    }

    private void loadCutoffSensitiveData(final int numberOfRatiosInChr1)
    {
        ratios.clear();

        loadChr1Data(numberOfRatiosInChr1);
        ratios.put(_2, ratio(_2, 1, 0.5));
        ratios.put(_2, ratio(_2, 1001, 0.48));
        ratios.put(_2, ratio(_2, 2001, 0.52));
        ratios.put(_2, ratio(_2, 3001, 0.51));
        ratios.put(_2, ratio(_2, 4001, 0.49));
        ratios.put(_2, ratio(_2, 5001, 0.50));
        ratios.put(_2, ratio(_2, 6001, 10.0));
        ratios.put(_2, ratio(_2, 10001, 10.0));
        ratios.put(_2, ratio(_2, 11001, 10.0));
        ratios.put(_2, ratio(_2, 12001, 10.7));
        ratios.put(_2, ratio(_2, 13001, 10.7));
        ratios.put(_2, ratio(_2, 14001, 0.01));
        ratios.put(_2, ratio(_2, 15001, 0.02));
        ratios.put(_2, ratio(_2, 16001, 0.01));
        ratios.put(_2, ratio(_2, 17001, 0.02));
    }

    private void loadChr1Data(final int ratioCount)
    {
        for(int i = 1; i < ratioCount; i++)
        {
            ratios.put(_1, ratio(_1, 1000 * i + 1, 0.5));
        }
    }

    private PositionDepth ratio(HumanChromosome chromosome, int start, double value)
    {
        return new PositionDepth(start, chromosome.shortName(), value);
    }
}
record PositionDepth(int position, String chromosome, double depth) implements GenomePosition
{
}
class PAS extends PerArmSegmenter<PositionDepth>
{
    protected PAS(final ListMultimap<Chromosome, PositionDepth> ratios,
            final ChrArmLocator chrArmLocator, final double gamma, final int threshold)
    {
        super(ratios, chrArmLocator, gamma, threshold);
    }

    @Override
    public double value(final PositionDepth ratio)
    {
        return ratio.depth();
    }

    @Override
    public DataForSegmentation buildSegmentationData(final List<PositionDepth> ratios)
    {
        double[] valuesForSegmentation = new double[ratios.size()];
        double[] rawValues = new double[ratios.size()];
        for(int i = 0; i < ratios.size(); i++)
        {
            PositionDepth ratio = ratios.get(i);
            rawValues[i] = ratio.depth();
            valuesForSegmentation[i] = ratio.depth();
        }
        return new DataForSegmentation(valuesForSegmentation, rawValues);
    }

    @Override
    public boolean isWindowed()
    {
        return true;
    }
}
