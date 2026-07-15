package com.hartwig.hmftools.common.segmentation;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._X;
import static com.hartwig.hmftools.common.segmentation.Arm.P;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.utils.pcf.PcfSegment;

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
    public void segmentationIsDeterministicTest() throws Exception
    {
        ratios.clear();
        int blockSize = 20;
        double shift = 300.0;
        for(HumanChromosome hc : HumanChromosome.values())
        {
            double multiplier = hc.ordinal() % 2 == 0 ? 1.0 : 100.0;
            double value = hc.intValue() * multiplier;
            int position = 1000;
            for(int i = 0; i < blockSize; i++)
            {
                double pointValue = value + (hc == _1 && i >= blockSize / 2 ? shift : 0.0);
                ratios.put(hc, ratio(hc, position, pointValue));
                position += 1000;
            }
        }

        PAS segmenter = new PAS(ratios, locator, GAMMA_100);
        Map<ChrArm, ChromosomeArmSegments<PositionDepth>> result0 = segmenter.getSegmentation(executor);

        // The key set iteration order previously affected the adjusted gamma and hence
        // the segmentation. We will run the algorithm with differently sorted key sets
        // and check the output compared to the first run.
        for(int i = 0; i < 100; i++)
        {
            List<Chromosome> toShuffle = new ArrayList<>(ratios.keySet().stream().toList());
            Collections.shuffle(toShuffle);
            ListMultimap<Chromosome, PositionDepth> copied = Multimaps.newListMultimap(new LinkedHashMap<>(), ArrayList::new);
            toShuffle.forEach(chromosome -> copied.putAll(chromosome, ratios.get(chromosome)));
            segmenter = new PAS(copied, locator, GAMMA_100);
            Map<ChrArm, ChromosomeArmSegments<PositionDepth>> result = segmenter.getSegmentation(executor);

            Assert.assertEquals(result0.keySet(), result.keySet());
            for(ChrArm chrArm : result.keySet())
            {
                ChromosomeArmSegments<PositionDepth> expected = result0.get(chrArm);
                ChromosomeArmSegments<PositionDepth> actual = result.get(chrArm);
                Assert.assertEquals(expected.Segments.size(), actual.Segments.size());
                for(int j = 0; j < expected.Segments.size(); j++)
                {
                    final PcfSegment expectedSegment = expected.Segments.get(j);
                    final PcfSegment actualSegment = actual.Segments.get(j);
                    Assert.assertEquals(expectedSegment.start(), actualSegment.start());
                    Assert.assertEquals(expectedSegment.end(), actualSegment.end());
                    Assert.assertEquals(expectedSegment.MeanRatio, actualSegment.MeanRatio, 0.0001);
                }
            }
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
            final ChrArmLocator chrArmLocator, final double gamma)
    {
        super(ratios, chrArmLocator, gamma);
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
