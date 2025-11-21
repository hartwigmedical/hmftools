package com.hartwig.hmftools.cobalt.segmentation;

import static com.hartwig.hmftools.cobalt.segmentation.Arm.P;
import static com.hartwig.hmftools.cobalt.segmentation.Arm.Q;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._X;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.nio.file.Files;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.pcf.CobaltSegment;
import com.hartwig.hmftools.common.utils.pcf.PCFFile;

import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class RatioSegmenterTest extends SegmentationTestBase
{
    private final ListMultimap<Chromosome, CobaltRatio> ratios = ArrayListMultimap.create();
    private ExecutorService executor;
    ChrArmLocator locator = cobaltRatio ->
    {
        Arm arm = cobaltRatio.position() > 20_000 ? Q : P;
        return new ChrArm(cobaltRatio.chr(), arm);
    };

    @Before
    public void setup()
    {
        executor = Executors.newFixedThreadPool(8);
        ratios.clear();
    }

    @After
    public void cleanup()
    {
        executor.shutdown();
    }

    @Test
    public void writeSegmentationFileTest() throws Exception
    {
        ratios.clear();
        // Since we're not testing the segmenter itself,
        // just use the same ratios for each chromosome.
        List<HumanChromosome> chromosomes = List.of(_1, _2, _3);
        chromosomes.forEach(chromosome ->
        {
            ratios.put(chromosome, ratioForValue(chromosome, 1, 0.5));
            ratios.put(chromosome, ratioForValue(chromosome, 1001, 0.51));
            ratios.put(chromosome, ratioForValue(chromosome, 2001, 0.49));
            ratios.put(chromosome, ratioForValue(chromosome, 3001, -1.0));
            ratios.put(chromosome, ratioForValue(chromosome, 4001, -1.0));
            ratios.put(chromosome, ratioForValue(chromosome, 5001, -1.0));
            ratios.put(chromosome, ratioForValue(chromosome, 6001, -1.0));
            ratios.put(chromosome, ratioForValue(chromosome, 13001, 10.09));
            ratios.put(chromosome, ratioForValue(chromosome, 14001, 9.92));
            ratios.put(chromosome, ratioForValue(chromosome, 15001, 10.02));
            ratios.put(chromosome, ratioForValue(chromosome, 16001, -1.0));

            ratios.put(chromosome, ratioForValue(chromosome, 20_001, 0.5));
            ratios.put(chromosome, ratioForValue(chromosome, 21_001, 0.51));
            ratios.put(chromosome, ratioForValue(chromosome, 22_001, 0.49));
            ratios.put(chromosome, ratioForValue(chromosome, 23_001, 0.52));
        });
        File tempDir = Files.createTempDirectory("rst").toFile();
        File outputFile = new File(tempDir, "rst.pcf");
        Assert.assertFalse(outputFile.exists());
        RatioSegmenter.writeTumorSegments(ratios, 100.0, V38, executor, outputFile.getAbsolutePath());

        ListMultimap<Chromosome, CobaltSegment> pcfData = PCFFile.readCobaltPcfFile(outputFile.getAbsolutePath());
        assertEquals(3, pcfData.keySet().size());
        List<CobaltSegment> regions1 = pcfData.get(_1);
        assertEquals(3, regions1.size());
        assertEquals(1, regions1.get(0).start());
        assertEquals(3000, regions1.get(0).end());
        double mean_1_0 = (tumorGCRatio(0) + tumorGCRatio(1) + tumorGCRatio(2)) / 3.0;
        assertEquals(mean_1_0, regions1.get(0).MeanRatio, 0.0001);
        assertEquals(13001, regions1.get(1).start());
        assertEquals(16000, regions1.get(1).end());
        double mean_1_1 = (tumorGCRatio(7) + tumorGCRatio(8) + tumorGCRatio(9)) / 3.0;
        assertEquals(mean_1_1, regions1.get(1).MeanRatio, 0.0001);
        assertEquals(20001, regions1.get(2).start());
        assertEquals(24000, regions1.get(2).end());
        double mean_1_2 = (tumorGCRatio(11) + tumorGCRatio(12) + tumorGCRatio(13) + tumorGCRatio(14)) / 4.0;
        assertEquals(mean_1_2, regions1.get(2).MeanRatio, 0.0001);
    }

    @Test
    public void handleMaskedWindowsWithinSegments() throws Exception
    {
        ratios.clear();

        ratios.put(_1, ratioForValue(_1, 1, 0.5));
        ratios.put(_1, ratioForValue(_1, 1001, -1.0));
        ratios.put(_1, ratioForValue(_1, 2001, -1.0));
        ratios.put(_1, ratioForValue(_1, 3001, 0.51));
        ratios.put(_1, ratioForValue(_1, 4001, 0.49));
        ratios.put(_1, ratioForValue(_1, 5001, -1.0));
        ratios.put(_1, ratioForValue(_1, 6001, -1.0));
        ratios.put(_1, ratioForValue(_1, 10001, 10.09));
        ratios.put(_1, ratioForValue(_1, 11001, -1.0));
        ratios.put(_1, ratioForValue(_1, 12001, -1.0));
        ratios.put(_1, ratioForValue(_1, 13001, 9.92));
        ratios.put(_1, ratioForValue(_1, 14001, -1.0));
        ratios.put(_1, ratioForValue(_1, 15001, 10.02));
        ratios.put(_1, ratioForValue(_1, 16001, -1.0));

        ratios.put(_1, ratioForValue(_1, 20_001, -1.0));
        ratios.put(_1, ratioForValue(_1, 21_001, 0.5));
        ratios.put(_1, ratioForValue(_1, 22_001, 0.51));
        ratios.put(_1, ratioForValue(_1, 23_001, 0.49));
        ratios.put(_1, ratioForValue(_1, 24_001, -1.0));
        ratios.put(_1, ratioForValue(_1, 25_001, -1.0));
        ratios.put(_1, ratioForValue(_1, 26_001, 0.52));
        File tempDir = Files.createTempDirectory("rst").toFile();
        File outputFile = new File(tempDir, "rst.pcf");
        Assert.assertFalse(outputFile.exists());
        RatioSegmenter.writeTumorSegments(ratios, 100.0, V38, executor, outputFile.getAbsolutePath());

        ListMultimap<Chromosome, CobaltSegment> pcfData = PCFFile.readCobaltPcfFile(outputFile.getAbsolutePath());
        assertEquals(1, pcfData.keySet().size());
        List<CobaltSegment> regions1 = pcfData.get(_1);
        assertEquals(3, regions1.size());
        assertEquals(1, regions1.get(0).start());
        assertEquals(5000, regions1.get(0).end());
        double mean_1_0 = (tumorGCRatio(0) + tumorGCRatio(3) + tumorGCRatio(4)) / 3.0;
        assertEquals(mean_1_0, regions1.get(0).MeanRatio, 0.0001);
        assertEquals(10001, regions1.get(1).start());
        assertEquals(16000, regions1.get(1).end());
        double mean_1_1 = (tumorGCRatio(7) + tumorGCRatio(10) + tumorGCRatio(12)) / 3.0;
        assertEquals(mean_1_1, regions1.get(1).MeanRatio, 0.0001);
        assertEquals(21001, regions1.get(2).start());
        assertEquals(27000, regions1.get(2).end());
        double mean_1_2 = (tumorGCRatio(15) + tumorGCRatio(16) + tumorGCRatio(17) + tumorGCRatio(20)) / 4.0;
        assertEquals(mean_1_2, regions1.get(2).MeanRatio, 0.0001);
    }

    @Test
    public void handleZeroes() throws Exception
    {
        ratios.clear();

        ratios.put(_1, ratio(_1, 1, 0.5));
        ratios.put(_1, ratio(_1, 1001, -1.0));
        ratios.put(_1, ratio(_1, 2001, -1.0));
        ratios.put(_1, ratio(_1, 3001, 0.51));
        ratios.put(_1, ratio(_1, 4001, 0.49));
        ratios.put(_1, ratio(_1, 5001, -1.0));
        ratios.put(_1, ratio(_1, 6001, -1.0));
        ratios.put(_1, ratio(_1, 10001, 100.09));
        ratios.put(_1, ratio(_1, 11001, 120.0));
        ratios.put(_1, ratio(_1, 12001, -1.0));
        ratios.put(_1, ratio(_1, 13001, 119.92));
        ratios.put(_1, ratio(_1, 14001, -1.0));
        ratios.put(_1, ratio(_1, 15001, 0.02));
        ratios.put(_1, ratio(_1, 16001, 0.0));
        ratios.put(_1, ratio(_1, 17001, 0.0));

        File tempDir = Files.createTempDirectory("rst").toFile();
        File outputFile = new File(tempDir, "rst.pcf");
        Assert.assertFalse(outputFile.exists());
        RatioSegmenter.writeTumorSegments(ratios, 100.0, V38, executor, outputFile.getAbsolutePath());

        ListMultimap<Chromosome, CobaltSegment> pcfData = PCFFile.readCobaltPcfFile(outputFile.getAbsolutePath());
        assertEquals(1, pcfData.keySet().size());
        List<CobaltSegment> regions1 = pcfData.get(_1);
        assertEquals(4, regions1.size());
        assertEquals(0.0, regions1.get(3).MeanRatio, 0.0001);
    }

    @Test
    public void chromosomeLabels() throws Exception
    {
        ratios.clear();
        List<HumanChromosome> chromosomes = List.of(_1, _2, _X);
        chromosomes.forEach(chromosome ->
        {
            ratios.put(chromosome, ratioForValue(chromosome, 1, 0.5));
            ratios.put(chromosome, ratioForValue(chromosome, 1001, 0.51));
            ratios.put(chromosome, ratioForValue(chromosome, 2001, 0.49));
        });
        File tempDir = Files.createTempDirectory("rst").toFile();
        File outputFile = new File(tempDir, "rst.38.pcf");
        Assert.assertFalse(outputFile.exists());
        RatioSegmenter.writeTumorSegments(ratios, 100.0, V38, executor, outputFile.getAbsolutePath());

        ListMultimap<Chromosome, CobaltSegment> pcfData = PCFFile.readCobaltPcfFile(outputFile.getAbsolutePath());
        assertEquals(3, pcfData.keySet().size());
        assertEquals(V38.versionedChromosome("1"), pcfData.get(_1).get(0).chromosome());
        assertEquals(V38.versionedChromosome("2"), pcfData.get(_2).get(0).chromosome());
        assertEquals(V38.versionedChromosome("X"), pcfData.get(_X).get(0).chromosome());

        File outputFile37 = new File(tempDir, "rst.37.pcf");
        Assert.assertFalse(outputFile37.exists());
        RatioSegmenter.writeTumorSegments(ratios, 100.0, V37, executor, outputFile.getAbsolutePath());

        ListMultimap<Chromosome, CobaltSegment> pcfData37 = PCFFile.readCobaltPcfFile(outputFile.getAbsolutePath());
        assertEquals(3, pcfData37.keySet().size());
        assertEquals(V37.versionedChromosome("1"), pcfData37.get(_1).get(0).chromosome());
        assertEquals(V37.versionedChromosome("2"), pcfData37.get(_2).get(0).chromosome());
        assertEquals(V37.versionedChromosome("X"), pcfData37.get(_X).get(0).chromosome());
    }

    @Test
    public void referenceSegmenterTest()
    {
        ratios.put(_1, ratioForValue(_1, 1001, 1.0));
        ratios.put(_1, ratioForValue(_1, 2001, 1.0));
        RatioSegmenter segmenter = new ReferenceRatioSegmenter(ratios, locator, 50.0);
        CobaltRatio ratio = new CobaltRatio("1", 1, 2, 2.1, 2.2, 2.3, 3, 3.1, 3.2);
        assertEquals(ratio.referenceGCDiploidRatio(), segmenter.value(ratio), 0.00001);
    }

    @Test
    public void tumorSegmenterTest()
    {
        ratios.put(_1, ratioForValue(_1, 1001, 1.0));
        ratios.put(_1, ratioForValue(_1, 2001, 1.0));
        RatioSegmenter segmenter = new TumorRatioSegmenter(ratios, locator, 50.0);
        CobaltRatio ratio = new CobaltRatio("1", 1, 2, 2.1, 2.2, 2.3, 3, 3.1, 3.2);
        assertEquals(ratio.tumorGCRatio(), segmenter.value(ratio), 0.00001);
    }

    @Test
    public void useUniformPenaltyAcrossAllArmsIfThereAreFewerThan100KRatios() throws Exception
    {
        ratios.clear();

        // These chr1 ratios together with the 15 for chr2 make 99_999 ratios.
        // This will result in a segmentation penalty of 1 applied to both
        // chromosomes. This will lead to 5 segments being created for chr2.
        for (int i = 1; i < 99_985; i++)
        {
            ratios.put(_1, ratio(_1, 1000 * i +1, 0.5));
        }
        ratios.put(_2, ratio(_2, 1, 0.5));
        ratios.put(_2, ratio(_2, 1001, 0.48));
        ratios.put(_2, ratio(_2, 2001, 0.52));
        ratios.put(_2, ratio(_2, 3001, 0.51));
        ratios.put(_2, ratio(_2, 4001, 0.49));
        ratios.put(_2, ratio(_2, 5001, 0.50));
        ratios.put(_2, ratio(_2, 6001, 100.0));
        ratios.put(_2, ratio(_2, 10001, 100.09));
        ratios.put(_2, ratio(_2, 11001, 120.0));
        ratios.put(_2, ratio(_2, 12001, 110.0));
        ratios.put(_2, ratio(_2, 13001, 119.92));
        ratios.put(_2, ratio(_2, 14001, 0.01));
        ratios.put(_2, ratio(_2, 15001, 0.02));
        ratios.put(_2, ratio(_2, 16001, 0.1));
        ratios.put(_2, ratio(_2, 17001, 0.4));

        File tempDir = Files.createTempDirectory("rst").toFile();
        File outputFile = new File(tempDir, "rst.pcf");
        Assert.assertFalse(outputFile.exists());
        RatioSegmenter.writeTumorSegments(ratios, 100.0, V38, executor, outputFile.getAbsolutePath());

        ListMultimap<Chromosome, CobaltSegment> pcfData = PCFFile.readCobaltPcfFile(outputFile.getAbsolutePath());
        assertEquals(2, pcfData.keySet().size());
        List<CobaltSegment> regions2 = pcfData.get(_2);
        assertEquals(5, regions2.size());
    }

    @Test
    public void usePerArmPenaltyIfThereAreMoreThan100KRatios() throws Exception
    {
        ratios.clear();

        // These chr1 ratios together with the 15 for chr2 make 100,000 ratios.
        // This will result in a separate penalty being calculated for each
        // chromosome arm. This means that chr2 becomes a single segment (note that the input gamma is 100).
        for (int i = 1; i < 99_986; i++)
        {
            ratios.put(_1, ratio(_1, 1000 * i +1, 0.5));
        }
        ratios.put(_2, ratio(_2, 1, 0.5));
        ratios.put(_2, ratio(_2, 1001, 0.48));
        ratios.put(_2, ratio(_2, 2001, 0.52));
        ratios.put(_2, ratio(_2, 3001, 0.51));
        ratios.put(_2, ratio(_2, 4001, 0.49));
        ratios.put(_2, ratio(_2, 5001, 0.50));
        ratios.put(_2, ratio(_2, 6001, 100.0));
        ratios.put(_2, ratio(_2, 10001, 100.09));
        ratios.put(_2, ratio(_2, 11001, 120.0));
        ratios.put(_2, ratio(_2, 12001, 110.0));
        ratios.put(_2, ratio(_2, 13001, 119.92));
        ratios.put(_2, ratio(_2, 14001, 0.01));
        ratios.put(_2, ratio(_2, 15001, 0.02));
        ratios.put(_2, ratio(_2, 16001, 0.1));
        ratios.put(_2, ratio(_2, 17001, 0.4));

        File tempDir = Files.createTempDirectory("rst").toFile();
        File outputFile = new File(tempDir, "rst.pcf");
        Assert.assertFalse(outputFile.exists());
        RatioSegmenter.writeTumorSegments(ratios, 100.0, V38, executor, outputFile.getAbsolutePath());

        ListMultimap<Chromosome, CobaltSegment> pcfData = PCFFile.readCobaltPcfFile(outputFile.getAbsolutePath());
        assertEquals(2, pcfData.keySet().size());
        List<CobaltSegment> regions2 = pcfData.get(_2);
        assertEquals(1, regions2.size());
    }

    private double tumorGCRatio(int position)
    {
        return ratios.get(HumanChromosome._1).get(position).tumorGCRatio();
    }
}
