package com.hartwig.hmftools.cobalt.e2e;

import static com.hartwig.hmftools.cobalt.CobaltConfig.PCF_GAMMA;
import static com.hartwig.hmftools.cobalt.CobaltConfig.TARGET_REGION_NORM_FILE;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
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

import com.hartwig.hmftools.cobalt.CobaltApplication;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.purple.Gender;

import org.apache.commons.io.FileUtils;
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

    /*
    Things to test:
     - sex chromosomes
     - GC normalisation
     - counting of reads that lie across windows
     - gc bucket smoothing
     - filtering of regions of extreme gc
     */

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

        // ratio -> ratio/relativeEnrichment: (1.0, 1.0, 1.0)/(2.0, 0.5, 1.0) = (0.5, 2.0, 1.0)
        // Normalised again by dividing by mean so that the mean is 1.0:
        // mean = 3.5/3 = 7/6, so (0.5, 2.0, 1.0_) => (1/7, 4/7, 2/7)*3.0
        // 0.5 / 7/6 = 3/7, 2.0 / 7/6 = 12/7, 1.0 / 7/6 = 6/7
        double oneSeventh = 1.0/7.0;
        assertEquals(3.0 * 1.0 * oneSeventh, ratios.get(1).tumorGCRatio(), 0.01);
        assertEquals(3.0 * 4.0 * oneSeventh, ratios.get(2).tumorGCRatio(), 0.01);
        assertEquals(3.0 * 2.0 * oneSeventh, ratios.get(3).tumorGCRatio(), 0.01);
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

    @Test
    public void readsAreApportionedToWindows() throws Exception
    {
        // 1 chr of length 6000
        // 1:1201-2201 depth 10
        // 1:2201-3201 depth 100
        // 1:3201-4201 depth 50
        sample = "three_windows_with_offset";
        bamFile = getBam(sample);
        regionOffset = 1_200;

        createStandardChr1GCFile(5_000);
        createStandardChr1PanelFile(5_000, 1.0001);

        runCobalt();

        // window 1001-2000 has 0.8 * 10 = 8
        // window 2001-3000 has 0.2 * 10 + 0.8 * 100 = 82
        // window 3001-4000 has 0.2 * 100 + 0.8 * 50 = 60
        // window 4001-5000 has 0.2 * 50 = 10
        List<CobaltRatio> ratios = ratioResults.get(_1);
        assertEquals(8.0, ratios.get(1).tumorReadDepth(), 0.01);
        assertEquals(82.0, ratios.get(2).tumorReadDepth(), 0.01);
        assertEquals(60.0, ratios.get(3).tumorReadDepth(), 0.01);
        assertEquals(10.0, ratios.get(4).tumorReadDepth(), 0.01);
    }

    @Test
    public void regionsOfExtremeGCAreFilteredOut() throws Exception
    {
        // 1 chr of length 101_000
        // 1:1001-2001 depth 1, gc ratio = 0.00
        // 1:2001-3001 depth 1, gc ratio = 0.01
        // 1:3001-4001 depth 1, gc ratio = 0.02
        // etc
        sample = "increasing_gc_per_window";
        bamFile = getBam(sample);
        regionOffset = 1000;

        createStandardChr1GCFile(101_000);
        createStandardChr1PanelFile(101_000, 1.0001);

        runCobalt();
        // Upper and lower limits for gc ratio are 0.26 and 0.68 respectively
        List<CobaltRatio> ratios = ratioResults.get(_1);
        for (int i=1; i < 26; i++)
        {
            assertEquals(-1.0, ratios.get(i).tumorGCRatio(), 0.01);
        }
        for (int i=26; i < 68; i++)
        {
            assertEquals(1.0, ratios.get(i).tumorGCRatio(), 0.01);
        }
        for (int i=69; i < ratios.size(); i++)
        {
            assertEquals(-1.0, ratios.get(i).tumorGCRatio(), 0.01);
        }
    }

    @Test
    public void filterOutWindowsWhereGcProfileDataIsMissing() throws Exception
    {
        // 1 chr of length 5000
        // 1:1001-4000 depth 100
        setupForThreeWindowBam();

        // Overwrite the gc profile
        gcProfile = new File(tempDir, "GC_profile.1000bp.38.cnp");
        GcProfilesUtilities gcFileWriter = new GcProfilesUtilities();
        gcFileWriter.addSection(new ConstantMappablePercentageGcFileSection(_1, regionOffset, regionOffset, 0.99));
        gcFileWriter.addSection(new ConstantMappablePercentageGcFileSection(_1, regionOffset + 2000, regionOffset + 2_000, 0.99));
        gcFileWriter.write(gcProfile);

        runCobalt();

        List<CobaltRatio> ratios = ratioResults.get(_1);
        assertEquals(1.0, ratios.get(1).tumorGCRatio(), 0.01);
        assertEquals(-1.0, ratios.get(2).tumorGCRatio(), 0.01);
        assertEquals(1.0, ratios.get(3).tumorGCRatio(), 0.01);
    }

    @Test
    public void gcWindowsAreSmoothed() throws Exception
    {
        // 1 chr of length 60_000
        // 1001-2000, 2001-3000, 3001-4000   : gc 0.40, depth 10
        // 4001-5000                         : gc 0.40, depth 20
        // 5001-6000                         : gc 0.40, depth 50
        // 6001-7000, 7001-8000, 8001-9000   : gc 0.41, depth 11
        // 9001-10_000                       : gc 0.41, depth 21
        // 10_001-11_000                     : gc 0.41, depth 51
        // ...
        // 51_001-52_000, 52_001-53_000, 53_001-54_000 : gc 0.50, depth 20
        // 54_001-55_000                               : gc 0.50, depth 30
        // 55_001-56_000                               : gc 0.50, depth 60
        sample = "gc_40_to_50";
        bamFile = getBam(sample);
        regionOffset = 1_000;

        createStandardChr1GCFile(60_000);
        createStandardChr1PanelFile(60_000, 1.000001);

        runCobalt();

        // Normalisation calculations
        // Step 1: get mean and median of non-zero depths and calculate median/mean.
        // There are 59 read windows. Apart from the 4 depth 0 windows (positions > 56000) the depths are:
        // 10, 10, 10, 11, 11, 11, ..., 19, 19, 19, 20, 20, 20, 20, 21, 22, ...,30, 50, 51, ..., 60
        // The median is the 28th, which is 19:
        // 1-3 10, 4-6 11, ... , 28-30 19
        // The total is  3 * (10 + 11 + ... + 20) + (20 + 21 + ... + 30) + (50 + 51 + ... + 60) = 3 * 165 + 275 + 605 = 1375
        // mean = 1375 / 55 = 25
        // GC median normalisation factor = median/mean = 19/25 = 0.76

        // Step 2: assign each window to a GC bucket. Calculate the bucket median depths and smooth these median values
        // gc buckets:
        // bucket median count
        // 40 10 5
        // 41 11 5
        // ...
        // 50 20 5
        // These get smoothed to:
        // 41 11 5
        // ...
        // 49 19 5
        // (so no change except the first and last values are removed).

        // Step 3: multiply each window's read depth by the median normalisation factor
        // and divide by the median depth for that gc value. Note that the median
        // depths are 0.41 11, 0.42 12,..., 0.49 19
        // 6001, 7001, 8001: 11 -> 11*0.76/11
        // 9001 21 -> 21*0.76/11
        // 10_001 51 -> 51*0.76/11
        // 11_001, 12_001, 13_001: 12 -> 12*0.76/12
        // 14_001 22 -> 22*0.76/12
        // 15_001 51 -> 53*0.76/12
        // ...

        // Step 4: get the median of the on-target windows and normalise by dividing by this.
        // The median is 0.76 so the
        // 6001, 7001, 8001: 11 -> 11/11
        // 9001 21 -> 21/11
        // 10_001 51 -> 51/11
        // 11_001, 12_001, 13_001: 12 -> 12/12
        // 14_001 22 -> 22/12
        // 15_001 51 -> 53/12
        // ...
        // ????????????????? We've now multiplied by the median and then divided by it.
        // Maybe we should not do either of these steps.

        // Step 5: normalise by dividing by the mean of the non-negative values
        // so that the mean of the result is 1.0.
        // The values are:
        // (11, 11, 11, 21, 51)/11, (12, 12, 12, 22, 52)/12, ..., (19, 19, 19, 29, 59)/19
        // Total = 105/11 + 110/12 + ... + 145/19 = 75.94
        // Count = 45, so mean = 1.688
        double mean = 1.688;

        List<CobaltRatio> ratios = ratioResults.get(_1);
        checkTumorRatio(-1.0, 0, 1, 2, 3, 4,5 );
        checkTumorRatio(1.0/mean, 6,7,8 );
    }

    private void checkTumorRatio(double expected, int...indices)
    {
        List<CobaltRatio> ratios = ratioResults.get(_1);
        for(final int index : indices)
        {
            assertEquals(expected, ratios.get(index).tumorGCRatio(), 0.01);
        }
    }

    private void setupForSingleWindowBam() throws IOException
    {
        sample = "one_window";
        bamFile = getBam(sample);
        regionOffset = 1_000;

        createStandardChr1GCFile(2_000);
        createStandardChr1PanelFile(2_000, 1.00);
    }

    private void setupForThreeWindowBam() throws IOException
    {
        // 1 chr of length 5000
        // 1:1001-4000 depth 100
        sample = "three_windows";
        bamFile = getBam(sample);
        regionOffset = 1_000;

        createStandardChr1GCFile(4_000);
        createStandardChr1PanelFile(4_000, 1.0000001);
    }

    private void createStandardChr1PanelFile(final int length, final double relativeEnrichment) throws IOException
    {
        panelNormalisation = new File(tempDir, "ThePanel.tsv");
        PanelFileWriter panelWriter = new PanelFileWriter();
        panelWriter.addSection(new PanelFileSection(_1, regionOffset, regionOffset + length, relativeEnrichment));
        panelWriter.write(panelNormalisation);
    }

    private void createStandardChr1GCFile(final int length) throws IOException
    {
        gcProfile = new File(tempDir, "GC_profile.1000bp.38.cnp");
        GcProfilesUtilities gcFileWriter = new GcProfilesUtilities();
        gcFileWriter.addSection(new ConstantGcFileSection(_1, regionOffset, regionOffset + length, 0.5));
        gcFileWriter.write(gcProfile);
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
