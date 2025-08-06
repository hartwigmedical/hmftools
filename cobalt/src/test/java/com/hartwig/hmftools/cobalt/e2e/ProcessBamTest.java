package com.hartwig.hmftools.cobalt.e2e;

import static com.hartwig.hmftools.cobalt.CobaltConfig.PCF_GAMMA;
import static com.hartwig.hmftools.cobalt.CobaltConfig.TARGET_REGION_NORM_FILE;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.gc.GCProfile.MIN_MAPPABLE_PERCENTAGE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.CobaltApplication;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.purple.Gender;

import org.apache.commons.io.FileUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class ProcessBamTest
{
    private File tempDir;
    private String sample;
    private int regionOffset;
    private File bamFile;
    private File gcProfile;
    private File panelNormalisation;
    private File outputDir;
    private Map<Chromosome, List<CobaltRatio>> ratioResults;

    @Before
    public void setup() throws Exception
    {
        tempDir = new File("/Users/timlavers/work/junk/rubbish"); // TODO

        outputDir = new File(tempDir, "output");
        outputDir.mkdirs();
        FileUtils.cleanDirectory(outputDir);
    }

    @Test
    public void singleWindow() throws Exception
    {
        setupForSingleWindowBam();

        runCobalt();

        assertEquals(1, ratioResults.size());
        List<CobaltRatio> ratios = ratioResults.get(_1);
        assertEquals(3, ratios.size());
        assertEquals(1, ratios.get(0).position());
        assertEquals(0.0, ratios.get(0).tumorReadDepth(), 0.01);
        assertEquals(-1.0, ratios.get(0).tumorGCRatio(), 0.01);

        assertEquals(1001, ratios.get(1).position());
        assertEquals(100.0, ratios.get(1).tumorReadDepth(), 0.01);
        assertEquals(1.0, ratios.get(1).tumorGCRatio(), 0.01);
        assertEquals(0.5, ratios.get(1).tumorGcContent(), 0.01);

        assertEquals(2001, ratios.get(2).position());
        assertEquals(0.0, ratios.get(2).tumorReadDepth(), 0.01);
        assertEquals(-1.0, ratios.get(2).tumorGCRatio(), 0.01);
    }

    @Test
    public void panelNormalisationIsApplied() throws Exception
    {
        // 1 chr of length 5000
        // 1:1001-4000 depth 100
        setupForThreeWindowBam();

        // Overwrite the normalisation profile
        panelNormalisation = new File(tempDir, "ThePanel.tsv");
        PanelFileWriter panelWriter = new PanelFileWriter();
        panelWriter.addSection(new PanelFileSection(_1, regionOffset, regionOffset, 2.000));
        panelWriter.addSection(new PanelFileSection(_1, regionOffset + 1000, regionOffset + 1000, 0.5000));
        panelWriter.addSection(new PanelFileSection(_1, regionOffset + 2000, regionOffset + 2000, 1.0000));
        panelWriter.write(panelNormalisation);

        runCobalt();
        List<CobaltRatio> ratios = ratioResults.get(_1);

        assertEquals(0.5, ratios.get(1).tumorGCRatio(), 0.01);
        assertEquals(2.0, ratios.get(2).tumorGCRatio(), 0.01);
        assertEquals(1.0, ratios.get(3).tumorGCRatio(), 0.01);
    }

    @Test
    public void gcProfileIsUsedToFilterOutNonMappableWindows() throws Exception
    {
        // 1 chr of length 5000
        // 1:1001-4000 depth 100
        setupForThreeWindowBam();

        // Overwrite the gc profile
        double mappability0 = MIN_MAPPABLE_PERCENTAGE + 0.01;
        double mappability2 = MIN_MAPPABLE_PERCENTAGE - 0.01;
        gcProfile = new File(tempDir, "GC_profile.1000bp.38.cnp");
        GcProfilesUtilities gcFileWriter = new GcProfilesUtilities();
        gcFileWriter.addSection(new ConstantMappablePercentageGcFileSection(_1, regionOffset, regionOffset, mappability0));
        gcFileWriter.addSection(new ConstantMappablePercentageGcFileSection(_1, regionOffset+ 1000, regionOffset + 1_000, MIN_MAPPABLE_PERCENTAGE));
        gcFileWriter.addSection(new ConstantMappablePercentageGcFileSection(_1, regionOffset+ 2000, regionOffset + 2_000, mappability2));
        gcFileWriter.write(gcProfile);

        runCobalt();

        List<CobaltRatio> ratios = ratioResults.get(_1);
        assertEquals(1.0, ratios.get(1).tumorGCRatio(), 0.01);
        assertEquals(1.0, ratios.get(2).tumorGCRatio(), 0.01);
        assertEquals(-1.0, ratios.get(3).tumorGCRatio(), 0.01);
    }

    private void setupForSingleWindowBam() throws IOException
    {
        // 1 chr of length 3000
        // 1:1001-2000 depth 100
        sample = "one_window";
        bamFile = getBam(sample);
        regionOffset = 1_000;

        gcProfile = new File(tempDir, "GC_profile.1000bp.38.cnp");
        GcProfilesUtilities gcFileWriter = new GcProfilesUtilities();
        gcFileWriter.addSection(new ConstantGcFileSection(_1, regionOffset, regionOffset + 2_000, 0.5));
        gcFileWriter.write(gcProfile);

        panelNormalisation = new File(tempDir, "ThePanel.tsv");
        PanelFileWriter panelWriter = new PanelFileWriter();
        panelWriter.addSection(new PanelFileSection(_1, regionOffset, regionOffset + 2_000, 1.00));
        panelWriter.write(panelNormalisation);
    }

    private void setupForThreeWindowBam() throws IOException
    {
        // 1 chr of length 5000
        // 1:1001-4000 depth 100
        sample = "three_windows";
        bamFile = getBam(sample);
        regionOffset = 1_000;

        gcProfile = new File(tempDir, "GC_profile.1000bp.38.cnp");
        GcProfilesUtilities gcFileWriter = new GcProfilesUtilities();
        gcFileWriter.addSection(new ConstantGcFileSection(_1, regionOffset, regionOffset + 4_000, 0.5));
        gcFileWriter.write(gcProfile);

        panelNormalisation = new File(tempDir, "ThePanel.tsv");
        PanelFileWriter panelWriter = new PanelFileWriter();
        panelWriter.addSection(new PanelFileSection(_1, regionOffset, regionOffset + 4_000, 1.00));
        panelWriter.write(panelNormalisation);
    }

    private File getBam(String sample)
    {
        return Path.of("src", "test", "resources", "bam", sample + ".sorted.bam").toFile();
    }

    private void runCobalt() throws Exception
    {
        String[] args = new String[12];
        int index = 0;
        args[index++] = String.format("-%s", TUMOR);
        args[index++] = String.format("%s", sample);
        args[index++] = String.format("-%s", TUMOR_BAM);
        args[index++] = String.format("%s", bamFile.getAbsolutePath());
        args[index++] = String.format("-%s", GC_PROFILE);
        args[index++] = String.format("%s", gcProfile.getAbsolutePath());
        args[index++] = String.format("-%s", TARGET_REGION_NORM_FILE);
        args[index++] = String.format("%s", panelNormalisation.getAbsolutePath());
        args[index++] = String.format("-%s", PCF_GAMMA);
        args[index++] = String.format("%d", 1);
        args[index++] = String.format("-%s", OUTPUT_DIR);
        args[index] = String.format("%s", outputDir.getAbsolutePath());

        CobaltApplication.main(args);

        File ratioFile = new File(outputDir, sample + ".cobalt.ratio.tsv.gz");
        assertTrue(ratioFile.exists());
        assertTrue(ratioFile.isFile());
        ratioResults = CobaltRatioFile.readWithGender(ratioFile.getAbsolutePath(), Gender.FEMALE, true);
    }

    @Test
    public void basic() throws Exception
    {
        // Example0 BAM has:
        // 1:10_000_001-10_100_000 depth 100
        // 1:10_200_001-10_300_000 depth 100
        // 1:10_400_001-10_500_000 depth 100
        // 2:10_000_001-10_100_000 depth 100
        // 2:10_200_001-10_300_000 depth 100
        // 2:10_400_001-10_500_000 depth 100

        // Example2 BAM has:
        // 1:10_000_001-10_100_000 depth 100
        // 1:10_200_001-10_300_000 depth 100
        // 1:10_400_001-10_500_000 depth 100
        // 2:10_000_001-10_100_000 depth 50
        // 2:10_200_001-10_300_000 depth 50
        // 2:10_400_001-10_500_000 depth 50

        // Example3 BAM has:
        // 1:10_000_001-10_200_000 depth 100
        // 1:10_200_001-10_400_000 depth 10
        // 1:10_400_001-10_600_000 depth 100

        // Example4 BAM has same form as Example3 but has random bases.

        // Example6
        // 1 chr of length 12000
        //1:2001-5000 depth 100
        //1:5001-8000 depth 10
        //1:8001-11000 depth 100

        // Example7
        // 1 chr of length 3000
        //1:1001-2000 depth 100
        String sample = "Example7";
        File tempDir = new File("/Users/timlavers/work/junk/rubbish");
        File bamFile = new File(tempDir, sample + ".sorted.bam");
        final int regionOffset = 1_000;

        RefGenomeInterface genome = new CobaltTestGenome();
        File gcProfile = new File(tempDir, "GC_profile.1000bp.38.cnp");
        GcProfilesUtilities gcFileWriter = new GcProfilesUtilities();
//        gcFileWriter.addSection(new GcFileSection(_1, regionOffset, regionOffset + 40_000, genome));
        gcFileWriter.addSection(new ConstantGcFileSection(_1, regionOffset, regionOffset + 2_000, 0.5));
//        gcFileWriter.addSection(new GcFileSection(_2, regionOffset, regionOffset + 600_000, genome));
        gcFileWriter.write(gcProfile);

        File panelNormalisation = new File(tempDir, "ThePanel.tsv");
        PanelFileWriter panelWriter = new PanelFileWriter();
        panelWriter.addSection(new PanelFileSection(_1, regionOffset, regionOffset + 2_000, 1.00));
//        panelWriter.addSection(new PanelFileSection(_1, regionOffset + 200_000, regionOffset + 300_000, 1.001));
//        panelWriter.addSection(new PanelFileSection(_1, regionOffset + 400_000, regionOffset + 500_000, 1.001));
//        panelWriter.addSection(new PanelFileSection(_2, regionOffset, regionOffset + 600_000, 1.001));
//        panelWriter.addSection(new PanelFileSection(_2, regionOffset + 200_000, regionOffset + 300_000, 1.001));
//        panelWriter.addSection(new PanelFileSection(_2, regionOffset + 400_000, regionOffset + 500_000, 1.001));
        panelWriter.write(panelNormalisation);

        File outputDir = new File(tempDir, "output");
        outputDir.mkdirs();
        FileUtils.cleanDirectory(outputDir);

        File debugDir = new File(tempDir, "debug");
        debugDir.mkdirs();
        FileUtils.cleanDirectory(debugDir);

        String[] args = new String[12];
        int index = 0;
        args[index++] = String.format("-%s", TUMOR);
        args[index++] = String.format("%s", sample);
        args[index++] = String.format("-%s", TUMOR_BAM);
        args[index++] = String.format("%s", bamFile.getAbsolutePath());
        args[index++] = String.format("-%s", GC_PROFILE);
        args[index++] = String.format("%s", gcProfile.getAbsolutePath());
        args[index++] = String.format("-%s", TARGET_REGION_NORM_FILE);
        args[index++] = String.format("%s", panelNormalisation.getAbsolutePath());
        args[index++] = String.format("-%s", "pcf_gamma");
        args[index++] = String.format("%d", 1);
        args[index++] = String.format("-%s", OUTPUT_DIR);
        args[index] = String.format("%s", outputDir.getAbsolutePath());

        CobaltApplication.main(args);

        System.out.println("Output directory: " + outputDir.getAbsolutePath());

        File[] filesWritten = outputDir.listFiles();
        System.out.println("Written " + filesWritten.length + " files");

        File ratioFile = new File(outputDir, sample + ".cobalt.ratio.tsv.gz");
        assertTrue(ratioFile.exists());
        assertTrue(ratioFile.isFile());
        Map<Chromosome, List<CobaltRatio>> ratioResults = CobaltRatioFile.readWithGender(ratioFile.getAbsolutePath(), Gender.FEMALE, true);

        assertEquals(1, ratioResults.size());
/*
chromosome	position	referenceReadDepth	tumorReadDepth	referenceGCRatio	tumorGCRatio	referenceGCDiploidRatio	referenceGCContent	tumorGCContent
chr1	1	-1	0	-1	-1	-1	-1	-1
chr1	1001	-1	100	-1	1	-1	-1	0.5
chr1	2001	-1	0	-1	-1	-1	-1	-1
 */
    }
}
