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

public class RatioSegmenterTest
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
        double mean_1_0 = (tumorGCRatio(_1, 0) + tumorGCRatio(_1, 1) + tumorGCRatio(_1, 2)) / 3.0;
        assertEquals(mean_1_0, regions1.get(0).MeanRatio, 0.0001);
        assertEquals(13001, regions1.get(1).start());
        assertEquals(16000, regions1.get(1).end());
        double mean_1_1 = (tumorGCRatio(_1, 7) + tumorGCRatio(_1, 8) + tumorGCRatio(_1, 9)) / 3.0;
        assertEquals(mean_1_1, regions1.get(1).MeanRatio, 0.0001);
        assertEquals(20001, regions1.get(2).start());
        assertEquals(24000, regions1.get(2).end());
        double mean_1_2 = (tumorGCRatio(_1, 11) + tumorGCRatio(_1, 12) + tumorGCRatio(_1, 13) + tumorGCRatio(_1, 14)) / 4.0;
        assertEquals(mean_1_2, regions1.get(2).MeanRatio, 0.0001);
    }

    @Test
    public void handleMaskedWindowsWithinSegments() throws Exception
    {
        ratios.clear();

        ratios.put(_1, cr(_1, 1, 0.5));
        ratios.put(_1, cr(_1, 1001, -1.0));
        ratios.put(_1, cr(_1, 2001, -1.0));
        ratios.put(_1, cr(_1, 3001, 0.51));
        ratios.put(_1, cr(_1, 4001, 0.49));
        ratios.put(_1, cr(_1, 5001, -1.0));
        ratios.put(_1, cr(_1, 6001, -1.0));
        ratios.put(_1, cr(_1, 10001, 10.09));
        ratios.put(_1, cr(_1, 11001, -1.0));
        ratios.put(_1, cr(_1, 12001, -1.0));
        ratios.put(_1, cr(_1, 13001, 9.92));
        ratios.put(_1, cr(_1, 14001, -1.0));
        ratios.put(_1, cr(_1, 15001, 10.02));
        ratios.put(_1, cr(_1, 16001, -1.0));

        ratios.put(_1, cr(_1, 20_001, -1.0));
        ratios.put(_1, cr(_1, 21_001, 0.5));
        ratios.put(_1, cr(_1, 22_001, 0.51));
        ratios.put(_1, cr(_1, 23_001, 0.49));
        ratios.put(_1, cr(_1, 24_001, -1.0));
        ratios.put(_1, cr(_1, 25_001, -1.0));
        ratios.put(_1, cr(_1, 26_001, 0.52));
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
        double mean_1_0 = (tumorGCRatio(_1, 0) + tumorGCRatio(_1, 3) + tumorGCRatio(_1, 4)) / 3.0;
        assertEquals(mean_1_0, regions1.get(0).MeanRatio, 0.0001);
        assertEquals(10001, regions1.get(1).start());
        assertEquals(16000, regions1.get(1).end());
        double mean_1_1 = (tumorGCRatio(_1, 7) + tumorGCRatio(_1, 10) + tumorGCRatio(_1, 12)) / 3.0;
        assertEquals(mean_1_1, regions1.get(1).MeanRatio, 0.0001);
        assertEquals(21001, regions1.get(2).start());
        assertEquals(27000, regions1.get(2).end());
        double mean_1_2 = (tumorGCRatio(_1, 15) + tumorGCRatio(_1, 16) + tumorGCRatio(_1, 17) + tumorGCRatio(_1, 20)) / 4.0;
        assertEquals(mean_1_2, regions1.get(2).MeanRatio, 0.0001);
    }

    @Test
    public void handleZeroes() throws Exception
    {
        ratios.clear();

        ratios.put(_1, crDirect(_1, 1, 0.5));
        ratios.put(_1, crDirect(_1, 1001, -1.0));
        ratios.put(_1, crDirect(_1, 2001, -1.0));
        ratios.put(_1, crDirect(_1, 3001, 0.51));
        ratios.put(_1, crDirect(_1, 4001, 0.49));
        ratios.put(_1, crDirect(_1, 5001, -1.0));
        ratios.put(_1, crDirect(_1, 6001, -1.0));
        ratios.put(_1, crDirect(_1, 10001, 100.09));
        ratios.put(_1, crDirect(_1, 11001, 120.0));
        ratios.put(_1, crDirect(_1, 12001, -1.0));
        ratios.put(_1, crDirect(_1, 13001, 119.92));
        ratios.put(_1, crDirect(_1, 14001, -1.0));
        ratios.put(_1, crDirect(_1, 15001, 0.02));
        ratios.put(_1, crDirect(_1, 16001, 0.0));
        ratios.put(_1, crDirect(_1, 17001, 0.0));

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
            ratios.put(chromosome, cr(chromosome, 1, 0.5));
            ratios.put(chromosome, cr(chromosome, 1001, 0.51));
            ratios.put(chromosome, cr(chromosome, 2001, 0.49));
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

    private double tumorGCRatio(Chromosome chromosome, int position)
    {
        return ratios.get(chromosome).get(position).tumorGCRatio();
    }

    private static CobaltRatio cr(HumanChromosome chromosome, int start, double value)
    {
        // The RatioSegmenter uses the log of the values.
        // So that we know how the segmentation works out, we will anti-log them first.
        double v = value >= 0 ? Math.pow(2, value) : -1.0;
        return new CobaltRatio(chromosome.shortName(), start, v, v, v, v, v, v, v);
    }

    private static CobaltRatio crDirect(HumanChromosome chromosome, int start, double v)
    {
        return new CobaltRatio(chromosome.shortName(), start, v, v, v, v, v, v, v);
    }
}
