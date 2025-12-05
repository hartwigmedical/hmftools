package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.nio.file.Files;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;
import com.hartwig.hmftools.common.utils.pcf.PCFFile;
import com.hartwig.hmftools.common.utils.pcf.PcfSegment;

import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class BAFSegmenterTest
{
    private final RefGenomeVersion refGenomeVersion = RefGenomeVersion.V38;
    private final String chr1 = refGenomeVersion.versionedChromosome(HumanChromosome._1);
    private final String chr2 = refGenomeVersion.versionedChromosome(HumanChromosome._2);
    private final AmberBAF baf1_0 = new AmberBAF(chr1, 100, 0.5, 100, 0.51, 50);
    private final AmberBAF baf1_1 = new AmberBAF(chr1, 200, 0.51, 200, 0.52, 51);
    private final AmberBAF baf1_2 = new AmberBAF(chr1, 300, 0.49, 300, 0.53, 52);
    private final AmberBAF baf2_0 = new AmberBAF(chr2, 100, 0.5, 100, 0.51, 50);
    private final AmberBAF baf2_1 = new AmberBAF(chr2, 200, 0.51, 200, 0.52, 51);
    private final AmberBAF baf2_2 = new AmberBAF(chr2, 300, 0.49, 300, 0.53, 52);

    ExecutorService executor;
    File tempDir;

    @Before
    public void setup() throws Exception
    {
        executor = Executors.newFixedThreadPool(8);
        tempDir = Files.createTempDirectory("rst").toFile();
        FileUtils.cleanDirectory(tempDir);
    }

    @After
    public void cleanup()
    {
        executor.shutdown();
    }

    @Test
    public void segmentationTest() throws Exception
    {
        File outputFile = new File(tempDir, "rst.pcf");
        Assert.assertFalse(outputFile.exists());

        List<AmberBAF> bafList = List.of(baf1_0, baf1_1, baf1_2, baf2_0, baf2_1, baf2_2);
        BAFSegmenter.writeSegments(bafList, refGenomeVersion, executor, outputFile.getAbsolutePath());

        ListMultimap<Chromosome, PcfSegment> pcfData = PCFFile.readPcfFile(outputFile.getAbsolutePath());
        assertEquals(2, pcfData.keySet().size());
        List<PcfSegment> regions1 = pcfData.get(_1);
        assertEquals(1, regions1.size());
        assertEquals(100, regions1.get(0).start());
        assertEquals(300, regions1.get(0).end());
        double mean_1_0 = (baf1_0.tumorModifiedBAF() + baf1_1.tumorModifiedBAF() + baf1_2.tumorModifiedBAF()) / 3.0;
        assertEquals(mean_1_0, regions1.get(0).MeanRatio, 0.0001);
    }

    @Test
    public void isWindowed()
    {
        ListMultimap<Chromosome, AmberBAF> bafData = ArrayListMultimap.create();
        bafData.put(_1, baf1_0);
        BAFSegmenter segmenter = new BAFSegmenter(bafData, ChrArmLocator.defaultLocator(refGenomeVersion));
        Assert.assertFalse(segmenter.isWindowed());
    }

    @Test
    public void handleEmptyInput() throws Exception
    {
        File outputFile = new File(tempDir, "rst.pcf");
        Assert.assertFalse(outputFile.exists());

        List<AmberBAF> bafList = List.of();
        BAFSegmenter.writeSegments(bafList, refGenomeVersion, executor, outputFile.getAbsolutePath());

        ListMultimap<Chromosome, PcfSegment> pcfData = PCFFile.readPcfFile(outputFile.getAbsolutePath());
        assertEquals(0, pcfData.keySet().size());
    }
}
