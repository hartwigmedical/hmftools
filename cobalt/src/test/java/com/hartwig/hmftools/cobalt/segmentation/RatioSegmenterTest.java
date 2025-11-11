package com.hartwig.hmftools.cobalt.segmentation;

import static com.hartwig.hmftools.cobalt.segmentation.Arm.P;
import static com.hartwig.hmftools.cobalt.segmentation.Arm.Q;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.google.common.base.Stopwatch;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.utils.RawCobaltRatio;
import com.hartwig.hmftools.cobalt.utils.RawCobaltRatioFile;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.segmentation.PiecewiseConstantFit;
import com.hartwig.hmftools.common.utils.pcf.PCFFile;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class RatioSegmenterTest
{
    private final ListMultimap<Chromosome, CobaltRatio> ratios = ArrayListMultimap.create();
    private ExecutorService executor;
    ChrArmLocator locator = cobaltRatio ->
    {
        Arm arm = cobaltRatio.position() > 20_000 ? Q : P;
        return new ChrArm(cobaltRatio.chr(), arm);
    };

    public RatioSegmenterTest()
    {
        // Since we're not testing the segmenter itself,
        // just use the same ratios for each chromosome.
        List<HumanChromosome> chromosomes = List.of(_1, _2, _3);
        chromosomes.forEach(chromosome ->
        {
            ratios.put(chromosome, cr(chromosome, 1, 0.5));
            ratios.put(chromosome, cr(chromosome, 1001, 0.51));
            ratios.put(chromosome, cr(chromosome, 2001, 0.49));
            ratios.put(chromosome, cr(chromosome, 3001, -1.0));
            ratios.put(chromosome, cr(chromosome, 4001, -1.0));
            ratios.put(chromosome, cr(chromosome, 5001, -1.0));
            ratios.put(chromosome, cr(chromosome, 6001, -1.0));
            ratios.put(chromosome, cr(chromosome, 13001, 10.09));
            ratios.put(chromosome, cr(chromosome, 14001, 9.92));
            ratios.put(chromosome, cr(chromosome, 15001, 10.02));
            ratios.put(chromosome, cr(chromosome, 16001, -1.0));

            ratios.put(chromosome, cr(chromosome, 20_001, 0.5));
            ratios.put(chromosome, cr(chromosome, 21_001, 0.51));
            ratios.put(chromosome, cr(chromosome, 22_001, 0.49));
            ratios.put(chromosome, cr(chromosome, 23_001, 0.52));
        });
    }

    @Before
    public void setup()
    {
        executor = Executors.newFixedThreadPool(12);
    }

    @After
    public void cleanup()
    {
        executor.shutdown();
    }

//    @Test
    public void coloTest() throws Exception
    {
        ratios.clear();
        File rawCobaltRatiosFile = new File("/Users/timlavers/work/junk/cobalt/colo829/COLO829v003T.cobalt.ratio.tsv.gz");
        RawCobaltRatioFile rawResultsFile = new RawCobaltRatioFile(rawCobaltRatiosFile.getAbsolutePath());
        List<RawCobaltRatio> rawRatios = rawResultsFile.read();

        rawRatios.forEach(rawCobaltRatio ->
        {
            CobaltRatio ratio = rawCobaltRatio.toCobaltRatio();
            ratios.put(ratio.chr(), ratio);
        });
        System.out.println("Read " + ratios.size() + " ratios.");
        File resultsDir = new File("/Users/timlavers/work/junk");
        Stopwatch stopwatch = Stopwatch.createUnstarted();
        DescriptiveStatistics statistics = new DescriptiveStatistics();
        for (int i=0; i<3; i++)
        {
            System.out.println("Run " + i);
            File outputFile = new File(resultsDir, "colo." + i + ".pcf");
            stopwatch.start();
            RatioSegmenter.writeTumorSegments(ratios, 100.0, V38, executor, outputFile.getAbsolutePath());
            stopwatch.stop();
            statistics.addValue(stopwatch.elapsed().toMillis());
            stopwatch.reset();
        }
        System.out.println("Average: " + statistics.getMean());
        System.out.println("Stats: " + statistics);
    }

    @Test
    public void writeSegmentationFileTest() throws Exception
    {
        File tempDir = Files.createTempDirectory("rst").toFile();
        File outputFile = new File(tempDir, "rst.pcf");
        Assert.assertFalse(outputFile.exists());
        RatioSegmenter.writeTumorSegments(ratios, 100.0, V38, executor, outputFile.getAbsolutePath());

        ListMultimap<Chromosome, ChrBaseRegion> pcfData = PCFFile.readCobaltPcfFile(outputFile.getAbsolutePath());
        assertEquals(3, pcfData.keySet().size());
    }

    @Test
    public void segmentationTest() throws Exception
    {
        RatioSegmenter segmenter = new TumorRatioSegmenter(ratios, locator, 50.0);
        Map<ChrArm, PiecewiseConstantFit> segments = segmenter.getSegmentation(executor);
        assertEquals(6, segments.size());
        final PiecewiseConstantFit pFit = new PiecewiseConstantFit(new int[] { 3, 3 }, new int[] { 0, 3 }, new double[] { 0.5, 10.01 });
        assertEquals(pFit, segments.get(ca(_1, P)));
        assertEquals(pFit, segments.get(ca(_2, P)));
        assertEquals(pFit, segments.get(ca(_3, P)));
        final PiecewiseConstantFit qFit = new PiecewiseConstantFit(new int[] { 4 }, new int[] { 0 }, new double[] { 0.505 });
        assertEquals(qFit, segments.get(ca(_1, Q)));
        assertEquals(qFit, segments.get(ca(_2, Q)));
        assertEquals(qFit, segments.get(ca(_3, Q)));
    }

    @Test
    public void referenceSegmenterTest()
    {
        RatioSegmenter segmenter = new ReferenceRatioSegmenter(ratios, locator, 50.0);
        CobaltRatio ratio = new CobaltRatio("1", 1, 2, 2.1, 2.2, 2.3, 3, 3.1, 3.2);
        assertEquals(ratio.referenceGCDiploidRatio(), segmenter.value(ratio), 0.00001);
    }

    @Test
    public void tumorSegmenterTest()
    {
        RatioSegmenter segmenter = new TumorRatioSegmenter(ratios, locator, 50.0);
        CobaltRatio ratio = new CobaltRatio("1", 1, 2, 2.1, 2.2, 2.3, 3, 3.1, 3.2);
        assertEquals(ratio.tumorGCRatio(), segmenter.value(ratio), 0.00001);
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
