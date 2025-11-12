package com.hartwig.hmftools.cobalt.e2e;

import static com.hartwig.hmftools.cobalt.CobaltConfig.PCF_GAMMA;
import static com.hartwig.hmftools.cobalt.CobaltConfig.TARGET_REGION_NORM_FILE;
import static com.hartwig.hmftools.cobalt.CobaltConfig.TUMOR_ONLY_DIPLOID_BED;
import static com.hartwig.hmftools.cobalt.CobaltConfig.USE_NEW_SEGMENTER;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._15;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._16;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._X;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._Y;
import static com.hartwig.hmftools.common.genome.gc.GCProfile.MIN_MAPPABLE_PERCENTAGE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.CobaltApplication;
import com.hartwig.hmftools.cobalt.testutils.DiploidRegionsFileWriter;
import com.hartwig.hmftools.cobalt.utils.RawCobaltRatio;
import com.hartwig.hmftools.cobalt.utils.RawCobaltRatioFile;
import com.hartwig.hmftools.common.cobalt.CobaltMedianRatioFile;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.pcf.CobaltSegment;
import com.hartwig.hmftools.common.utils.pcf.PCFFile;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;

import org.apache.commons.io.FileUtils;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class ProcessBamTest
{
    private File tempDir;
    private String sample;
    private String referenceSample;
    private int regionOffset;
    private File tumorBamFile;
    private File referenceBamFile;
    private File gcProfile;
    private File panelNormalisation;
    private File diploidBedFile;
    private File outputDir;
    private Map<Chromosome, List<CobaltRatio>> tumorRatioResults;
    private List<MedianRatio> medianRatioResults;
    private List<MedianRatio> referenceMedianRatioResults;

    @Before
    public void setup() throws Exception
    {
        sample = null;
        referenceSample = null;
        regionOffset = 0;
        tumorBamFile = null;
        referenceBamFile = null;
        tumorRatioResults = null;
        medianRatioResults = null;
        referenceMedianRatioResults = null;
        tempDir = Files.createTempDirectory("pbt").toFile();
        outputDir = new File(tempDir, "output");
        diploidBedFile = null;
        //noinspection ResultOfMethodCallIgnored
        outputDir.mkdirs();
        FileUtils.cleanDirectory(outputDir);
    }

    @Test
    public void tumorOnly() throws Exception
    {
        setupForSingleWindowBamTumorOnly();

        runCobalt();

        assertEquals(1, tumorRatioResults.size());
        List<CobaltRatio> ratios = tumorRatioResults.get(_1);
        assertEquals(3, ratios.size());
        assertEquals("chr1", ratios.get(0).chromosome());
        assertEquals(1, ratios.get(0).position());
        assertEquals(-1.0, ratios.get(0).referenceReadDepth(), 0.01);
        assertEquals(-1.0, ratios.get(0).referenceGCRatio(), 0.01);
        assertEquals(-1.0, ratios.get(0).referenceGCDiploidRatio(), 0.01);
        assertEquals(-1.0, ratios.get(0).referenceGcContent(), 0.01);
        assertEquals(0.0, ratios.get(0).tumorReadDepth(), 0.01);
        assertEquals(0.0, ratios.get(0).tumorGCRatio(), 0.01);
        assertEquals(0.5, ratios.get(0).tumorGcContent(), 0.01);

        assertEquals("chr1", ratios.get(1).chromosome());
        assertEquals(1001, ratios.get(1).position());
        assertEquals(-1.0, ratios.get(1).referenceReadDepth(), 0.01);
        assertEquals(-1.0, ratios.get(1).referenceGCRatio(), 0.01);
        assertEquals(-1.0, ratios.get(1).referenceGCDiploidRatio(), 0.01);
        assertEquals(-1.0, ratios.get(1).referenceGcContent(), 0.01);
        assertEquals(100.0, ratios.get(1).tumorReadDepth(), 0.01);
        assertEquals(1.0, ratios.get(1).tumorGCRatio(), 0.01);
        assertEquals(0.5, ratios.get(1).tumorGcContent(), 0.01);

        assertEquals(2001, ratios.get(2).position());
        assertEquals(-1.0, ratios.get(2).referenceReadDepth(), 0.01);
        assertEquals(-1.0, ratios.get(2).referenceGCRatio(), 0.01);
        assertEquals(-1.0, ratios.get(2).referenceGCDiploidRatio(), 0.01);
        assertEquals(-1.0, ratios.get(2).referenceGcContent(), 0.01);
        assertEquals(0.0, ratios.get(2).tumorReadDepth(), 0.01);
        assertEquals(0.0, ratios.get(2).tumorGCRatio(), 0.01);
        assertEquals(0.5, ratios.get(2).tumorGcContent(), 0.01);
    }

    @Test
    public void applyDiploidRegionFiltering() throws Exception
    {
        // chr1 and chr2 both of length 10_000
        // Each chr has each window of depth 100.
        sample = "chrs_1_and_2_10_windows_each";
        tumorBamFile = getBam(sample);
        regionOffset = 0;
        createStandardMultiChromosomeGCFile(10_000, _1, _2);
        createStandardMultipleChromosomePanelFile(10_000, 1.0001, _1, _2);

        // Create a diploid regions file that masks out the first
        // half of chr1 and the second half of chr2.
        diploidBedFile = new File(tempDir, "diploid_bed.tsv.gz");
        DiploidRegionsFileWriter diploidRegionsFileWriter = new DiploidRegionsFileWriter();
        diploidRegionsFileWriter.addSection(new DiploidFileSection(_1, 5000, 10000));
        diploidRegionsFileWriter.addSection(new DiploidFileSection(_2, 0, 5000));
        diploidRegionsFileWriter.write(diploidBedFile);
        runCobalt();

        assertEquals(2, tumorRatioResults.size());
        List<CobaltRatio> ratios1 = tumorRatioResults.get(_1);
        List<CobaltRatio> ratios2 = tumorRatioResults.get(_2);
        assertEquals(10, ratios1.size());
        for(int i = 0; i < 5; i++)
        {
            assertEquals(-1.0, ratios1.get(i).tumorGCRatio(), 0.01);
            assertEquals(1.0, ratios2.get(i).tumorGCRatio(), 0.01);
        }
        for(int i = 5; i < 10; i++)
        {
            assertEquals(1.0, ratios1.get(i).tumorGCRatio(), 0.01);
            assertEquals(-1.0, ratios2.get(i).tumorGCRatio(), 0.01);
        }
    }

    @Test
    public void germlineOnly() throws Exception
    {
        setupForSingleWindowBamGermlineOnly();

        runCobalt();

        assertEquals(1, tumorRatioResults.size());
        List<CobaltRatio> ratios = tumorRatioResults.get(_1);
        assertEquals(3, ratios.size());
        assertEquals("chr1", ratios.get(0).chromosome());
        assertEquals(1, ratios.get(0).position());
        assertEquals(0.0, ratios.get(0).referenceReadDepth(), 0.01);
        assertEquals(0.0, ratios.get(0).referenceGCDiploidRatio(), 0.01);
        assertEquals(0.0, ratios.get(0).referenceGCRatio(), 0.01);
        assertEquals(0.5, ratios.get(0).referenceGcContent(), 0.01);
        assertEquals(-1.0, ratios.get(0).tumorGCRatio(), 0.01);
        assertEquals(-1.0, ratios.get(0).tumorReadDepth(), 0.01);
        assertEquals(-1.0, ratios.get(0).tumorGcContent(), 0.01);

        assertEquals("chr1", ratios.get(1).chromosome());
        assertEquals(1001, ratios.get(1).position());
        assertEquals(100.0, ratios.get(1).referenceReadDepth(), 0.01);
        assertEquals(1.0, ratios.get(1).referenceGCRatio(), 0.01);
        assertEquals(1.0, ratios.get(1).referenceGCDiploidRatio(), 0.01);
        assertEquals(0.5, ratios.get(1).referenceGcContent(), 0.01);
        assertEquals(-1.0, ratios.get(1).tumorReadDepth(), 0.01);
        assertEquals(-1.0, ratios.get(1).tumorGCRatio(), 0.01);
        assertEquals(-1.0, ratios.get(1).tumorGcContent(), 0.01);

        assertEquals(2001, ratios.get(2).position());
        assertEquals(0.0, ratios.get(2).referenceReadDepth(), 0.01);
        assertEquals(0.0, ratios.get(2).referenceGCRatio(), 0.01);
        assertEquals(0.0, ratios.get(2).referenceGCDiploidRatio(), 0.01);
        assertEquals(0.5, ratios.get(2).referenceGcContent(), 0.01);
        assertEquals(-1.0, ratios.get(2).tumorReadDepth(), 0.01);
        assertEquals(-1.0, ratios.get(2).tumorGCRatio(), 0.01);
        assertEquals(-1.0, ratios.get(2).tumorGcContent(), 0.01);

        // Check the median ratios file.
        assertNull(medianRatioResults);
        double expectedMedianRatio = 100.0 / 33.333; // Smoothed median gc for bucket is 100/3

        assertEquals(1, referenceMedianRatioResults.size());
        assertEquals("chr1", referenceMedianRatioResults.get(0).Chromosome);
        assertEquals(expectedMedianRatio, referenceMedianRatioResults.get(0).MedianRatio, 0.0001);
        assertEquals(1, referenceMedianRatioResults.get(0).Count);
    }

    @Test
    public void tumorAndGermline() throws Exception
    {
        setupForSingleWindowBamTumorAndGermline();

        runCobalt();

        assertEquals(1, tumorRatioResults.size());
        List<CobaltRatio> ratios = tumorRatioResults.get(_1);
        assertEquals(3, ratios.size());
        assertEquals("chr1", ratios.get(0).chromosome());
        assertEquals(1, ratios.get(0).position());
        assertEquals(0.0, ratios.get(0).referenceReadDepth(), 0.01);
        assertEquals(0.0, ratios.get(0).referenceGCRatio(), 0.01);
        assertEquals(0.0, ratios.get(0).referenceGCDiploidRatio(), 0.01);
        assertEquals(0.5, ratios.get(0).referenceGcContent(), 0.01);
        assertEquals(0.0, ratios.get(0).tumorReadDepth(), 0.01);
        assertEquals(0.0, ratios.get(0).tumorGCRatio(), 0.01);
        assertEquals(0.5, ratios.get(0).tumorGcContent(), 0.01);

        assertEquals("chr1", ratios.get(1).chromosome());
        assertEquals(1001, ratios.get(1).position());
        assertEquals(100.0, ratios.get(1).referenceReadDepth(), 0.01);
        assertEquals(1.0, ratios.get(1).referenceGCRatio(), 0.01);
        assertEquals(1.0, ratios.get(1).referenceGCDiploidRatio(), 0.01);
        assertEquals(0.5, ratios.get(1).referenceGcContent(), 0.01);
        assertEquals(100.0, ratios.get(1).tumorReadDepth(), 0.01);
        assertEquals(1.0, ratios.get(1).tumorGCRatio(), 0.01);
        assertEquals(0.5, ratios.get(1).tumorGcContent(), 0.01);

        assertEquals(2001, ratios.get(2).position());
        assertEquals(0.0, ratios.get(2).referenceReadDepth(), 0.01);
        assertEquals(0.0, ratios.get(2).referenceGCRatio(), 0.01);
        assertEquals(0.0, ratios.get(2).referenceGCDiploidRatio(), 0.01);
        assertEquals(0.5, ratios.get(2).referenceGcContent(), 0.01);
        assertEquals(0.0, ratios.get(2).tumorReadDepth(), 0.01);
        assertEquals(0.0, ratios.get(2).tumorGCRatio(), 0.01);
        assertEquals(0.5, ratios.get(2).tumorGcContent(), 0.01);

        // Check the median ratios file.
        assertNull(medianRatioResults);
        double expectedMedianRatio = 100.0 / 33.333; // Smoothed median gc for bucket is 100/3

        assertEquals(1, referenceMedianRatioResults.size());
        assertEquals("chr1", referenceMedianRatioResults.get(0).Chromosome);
        assertEquals(expectedMedianRatio, referenceMedianRatioResults.get(0).MedianRatio, 0.0001);
        assertEquals(1, referenceMedianRatioResults.get(0).Count);
    }

    @Test
    public void useSuppliedGcIfCoverageIsLow() throws Exception
    {
        // 1 chr of length 4_000
        // 1:1000 depth 0.8 gc 40
        // 1001:2000 depth 0.8 gc 42
        // 2001:3000 depth 10 gc 44
        // 3001:4000 depth 10 gc 46
        // 4001:5000 depth 10 gc 48
        sample = "low_depth_regions";
        tumorBamFile = getBam(sample);
        regionOffset = 0;

        // Provide a GC profile that has GC 0.6 in the first two windows.
        gcProfile = new File(tempDir, "GC_profile.1000bp.38.cnp");
        GcProfilesUtilities gcFileWriter = new GcProfilesUtilities();
        gcFileWriter.addSection(new ConstantGcFileSection(_1, 1, 1001, 0.6));
        gcFileWriter.addSection(new ConstantGcFileSection(_1, 1001, 2001, 0.6));
        gcFileWriter.addSection(new ConstantGcFileSection(_1, 2001, 3001, 0.52));
        gcFileWriter.addSection(new ConstantGcFileSection(_1, 3001, 4001, 0.52));
        gcFileWriter.addSection(new ConstantGcFileSection(_1, 4001, 5001, 0.52));
        gcFileWriter.write(gcProfile);

        runCobalt(false);
        List<CobaltRatio> ratios = tumorRatioResults.get(_1);
        assertEquals(5, ratios.size());
        assertEquals(0.6, ratios.get(0).tumorGcContent(), 0.01); // depth 0
        assertEquals(0.6, ratios.get(1).tumorGcContent(), 0.01); // depth 0
        assertEquals(0.44, ratios.get(2).tumorGcContent(), 0.01); // depth 0
        assertEquals(0.46, ratios.get(3).tumorGcContent(), 0.01); // depth 0
        assertEquals(0.48, ratios.get(4).tumorGcContent(), 0.01); // depth 0
    }

    @Test
    public void panelNormalisationIsApplied() throws Exception
    {
        // 1 chr of length 5000
        // 1:1001-4000 depth 100
        setupForThreeWindowBamTumorOnly();

        // Overwrite the normalisation profile
        panelNormalisation = new File(tempDir, "ThePanel.tsv");
        PanelFileWriter panelWriter = new PanelFileWriter();
        panelWriter.addSection(new PanelFileSection(_1, 0, 0, 0.000));
        panelWriter.addSection(new PanelFileSection(_1, 1000, 1000, 2.000));
        panelWriter.addSection(new PanelFileSection(_1, 2000, 2000, 0.5000));
        panelWriter.addSection(new PanelFileSection(_1, 3000, 3000, 1.0000));
        panelWriter.addSection(new PanelFileSection(_1, 4000, 4000, 2.0000));
        panelWriter.write(panelNormalisation);

        runCobalt();
        List<CobaltRatio> ratios = tumorRatioResults.get(_1);

        // ratio -> ratio/relativeEnrichment: (1.0, 1.0, 1.0)/(2.0, 0.5, 1.0) = (0.5, 2.0, 1.0)
        // Normalised again by dividing by mean so that the mean is 1.0:
        // mean = 3.5/3 = 7/6, so (0.5, 2.0, 1.0_) => (1/7, 4/7, 2/7)*3.0
        // 0.5 / 7/6 = 3/7, 2.0 / 7/6 = 12/7, 1.0 / 7/6 = 6/7
        double oneSeventh = 1.0 / 7.0;
        assertEquals(3.0 * 1.0 * oneSeventh, ratios.get(1).tumorGCRatio(), 0.01);
        assertEquals(3.0 * 4.0 * oneSeventh, ratios.get(2).tumorGCRatio(), 0.01);
        assertEquals(3.0 * 2.0 * oneSeventh, ratios.get(3).tumorGCRatio(), 0.01);
    }

    @Test
    public void lowCoverage() throws Exception
    {
        // chr1 and chr2 each have length 100_000.
        // chr1 has depth 2 for each window, chr2 has depth 4.
        sample = "chr1_depth_2_chr2_depth_4";
        tumorBamFile = getBam(sample);
        regionOffset = 0;
        createStandardMultiChromosomeGCFile(100_000, _1, _2);
        runCobalt(false);

        assertEquals(2, tumorRatioResults.size());
        List<CobaltRatio> ratios1 = tumorRatioResults.get(_1);
        assertEquals(100, ratios1.size());
        List<CobaltRatio> ratios2 = tumorRatioResults.get(_2);
        assertEquals(100, ratios2.size());
        // The median depth is 3. Because this is <8, window consolidation is applied.
        // The number of windows merged into a single window is 80/3 rounded up to 30.
        // The middle of each consolidated list of windows is chosen as the representative position.
        // At these consolidation windows, the ratio is 2/3 or 4/3, depending on the chromosome.
        // At other positions, the ratio is -1.
        // Because there is no window for 105_000, a window is set at 95
        Set<Integer> representativeWindows = new HashSet<>();
        representativeWindows.add(14);
        representativeWindows.add(44);
        representativeWindows.add(74);
        representativeWindows.add(95);

        for(int i = 0; i < 100; i++)
        {
            if(representativeWindows.contains(i))
            {
                assertEquals(2.0, ratios1.get(i).tumorGCRatio(), 0.01);
                assertEquals(4.0, ratios2.get(i).tumorGCRatio(), 0.01);
            }
            else
            {
                assertEquals(-1.0, ratios1.get(i).tumorGCRatio(), 0.01);
                assertEquals(-1.0, ratios2.get(i).tumorGCRatio(), 0.01);
            }
        }
    }

    @Test
    public void lowCoverageTargeted() throws Exception
    {
        // chr1 and chr2 each have length 100_000.
        // chr1 has depth 2 for each window, chr2 has depth 4.
        sample = "chr1_depth_2_chr2_depth_4";
        tumorBamFile = getBam(sample);
        regionOffset = 0;
        createStandardMultiChromosomeGCFile(100_000, _1, _2);
        panelNormalisation = new File(tempDir, "ThePanel.tsv");
        // Apply panel normalisation that adjusts the low read counts in chr1.
        PanelFileWriter panelWriter = new PanelFileWriter();
        panelWriter.addSection(new PanelFileSection(_1, 0, 100_000, 0.501));
        panelWriter.addSection(new PanelFileSection(_2, 0, 100_000, 1.001));
        panelWriter.write(panelNormalisation);
        runCobalt();

        assertEquals(2, tumorRatioResults.size());
        List<CobaltRatio> ratios1 = tumorRatioResults.get(_1);
        assertEquals(100, ratios1.size());
        List<CobaltRatio> ratios2 = tumorRatioResults.get(_2);
        assertEquals(100, ratios2.size());
        for(int i = 0; i < 100; i++)
        {
            assertEquals(1.0, ratios1.get(i).tumorGCRatio(), 0.01);
            assertEquals(1.0, ratios2.get(i).tumorGCRatio(), 0.01);
        }
    }

    @Test
    public void borderlineMappableWindows() throws Exception
    {
        // 1 chr of length 5000
        // 1:1001-4000 depth 100
        setupForThreeWindowBamTumorOnly();

        // Overwrite the gc profile
        double mappability0 = MIN_MAPPABLE_PERCENTAGE + 0.01;
        double mappability2 = MIN_MAPPABLE_PERCENTAGE - 0.01;
        gcProfile = new File(tempDir, "GC_profile.1000bp.38.cnp");
        GcProfilesUtilities gcFileWriter = new GcProfilesUtilities();
        gcFileWriter.addSection(new ConstantMappablePercentageGcFileSection(_1, 0, 1000, mappability0));
        gcFileWriter.addSection(new ConstantMappablePercentageGcFileSection(_1, 2000, 2000, MIN_MAPPABLE_PERCENTAGE));
        gcFileWriter.addSection(new ConstantMappablePercentageGcFileSection(_1, 3000, 4_000, mappability2));
        gcFileWriter.write(gcProfile);

        runCobalt();

        List<CobaltRatio> ratios = tumorRatioResults.get(_1);
        assertEquals(0.0, ratios.get(0).tumorGCRatio(), 0.01); // depth 0
        assertEquals(1.0, ratios.get(1).tumorGCRatio(), 0.01);
        assertEquals(1.0, ratios.get(2).tumorGCRatio(), 0.01);
        assertEquals(-1.0, ratios.get(3).tumorGCRatio(), 0.01);
        assertEquals(-1.0, ratios.get(4).tumorGCRatio(), 0.01); // depth 0 but not mappable
    }

    @Test
    public void readsAreApportionedToWindows() throws Exception
    {
        // 1 chr of length 6000
        // 1:1201-2201 depth 10
        // 1:2201-3201 depth 100
        // 1:3201-4201 depth 50
        sample = "three_windows_with_offset";
        tumorBamFile = getBam(sample);

        createStandardChr1GCFile(6_000);
        createStandardChr1PanelFile(6_000, 1.0001);

        runCobalt();

        // window 1001-2000 has 0.8 * 10 = 8
        // window 2001-3000 has 0.2 * 10 + 0.8 * 100 = 82
        // window 3001-4000 has 0.2 * 100 + 0.8 * 50 = 60
        // window 4001-5000 has 0.2 * 50 = 10
        List<CobaltRatio> ratios = tumorRatioResults.get(_1);
        assertEquals(8.0, ratios.get(1).tumorReadDepth(), 0.01);
        assertEquals(82.0, ratios.get(2).tumorReadDepth(), 0.01);
        assertEquals(60.0, ratios.get(3).tumorReadDepth(), 0.01);
        assertEquals(10.0, ratios.get(4).tumorReadDepth(), 0.01);
    }

    @Test
    public void regionsOfExtremeGCAreFilteredOutInTargetedMode() throws Exception
    {
        // 1 chr of length 101_000
        // 1:1001-2001 depth 1, gc ratio = 0.00
        // 1:2001-3001 depth 1, gc ratio = 0.01
        // 1:3001-4001 depth 1, gc ratio = 0.02
        // etc
        sample = "increasing_gc_per_window_depth_1";
        tumorBamFile = getBam(sample);

        createStandardChr1GCFile(101_000);
        createStandardChr1PanelFile(101_000, 1.0001);

        runCobalt();
        // Upper and lower limits for gc ratio are 0.24 and 0.68 respectively
        List<CobaltRatio> ratios = tumorRatioResults.get(_1);
        for(int i = 1; i < 25; i++)
        {
            final CobaltRatio cobaltRatio = ratios.get(i);
            assertEquals(1000 * i + 1, ratios.get(i).position());
            assertEquals(0.01 * (i), ratios.get(i).tumorGcContent(), 0.01);
            assertEquals(-1.0, cobaltRatio.tumorGCRatio(), 0.01);
        }
        for(int i = 25; i < 68; i++)
        {
            assertEquals(1.0, ratios.get(i).tumorGCRatio(), 0.01);
        }
        for(int i = 68; i < ratios.size(); i++)
        {
            assertEquals(-1.0, ratios.get(i).tumorGCRatio(), 0.01);
        }
    }

    @Test
    public void filterOutUnmappableWindows() throws Exception
    {
        // 1 chr of length 5000
        // 1:1001-4000 depth 100
        setupForThreeWindowBamTumorOnly();

        // Overwrite the gc profile
        gcProfile = new File(tempDir, "GC_profile.1000bp.38.cnp");
        GcProfilesUtilities gcFileWriter = new GcProfilesUtilities();
        gcFileWriter.addSection(new ConstantMappablePercentageGcFileSection(_1, regionOffset, regionOffset, 0.99));
        gcFileWriter.addSection(new ConstantMappablePercentageGcFileSection(_1, regionOffset + 1000, regionOffset + 1000, 0.05));
        gcFileWriter.addSection(new ConstantMappablePercentageGcFileSection(_1, regionOffset + 2000, regionOffset + 2_000, 0.99));
        gcFileWriter.addSection(new ConstantMappablePercentageGcFileSection(_1, regionOffset + 3000, regionOffset + 3_000, 0.99));
        gcFileWriter.addSection(new ConstantMappablePercentageGcFileSection(_1, regionOffset + 4000, regionOffset + 4_000, 0.99));
        gcFileWriter.write(gcProfile);

        runCobalt();

        List<CobaltRatio> ratios = tumorRatioResults.get(_1);
        assertEquals(0.0, ratios.get(0).tumorGCRatio(), 0.01);
        assertEquals(-1.0, ratios.get(1).tumorGCRatio(), 0.01);
        assertEquals(1.0, ratios.get(2).tumorGCRatio(), 0.01);
        assertEquals(1.0, ratios.get(3).tumorGCRatio(), 0.01);
        assertEquals(0.0, ratios.get(4).tumorGCRatio(), 0.01);
    }

    @Test
    public void handleDepth0WindowsTargeted() throws Exception
    {
        // 1 chr of length 8000
        // 1:2000 depth 100, gc 0.5
        // 2001:4000 depth 100, gc 0.51
        // 4001:5000 depth 100, gc 0.52
        // 5001:6000 depth 0, gc 0
        // 6001:8000 depth 100, gc 0.53
        sample = "depth_0_windows";
        tumorBamFile = getBam(sample);
        regionOffset = 0;

        createStandardChr1GCFile(8000);
        createStandardChr1PanelFile(8000, 1.000001);

        runCobalt();

        List<CobaltRatio> ratios = tumorRatioResults.get(_1);
        assertEquals(4001, ratios.get(4).position());
        assertTrue(ratios.get(4).tumorGCRatio() > 0);
        assertEquals(100.0, ratios.get(4).tumorReadDepth(), 0.01);
        assertEquals(0.52, ratios.get(4).tumorGcContent(), 0.01);

        assertEquals(5001, ratios.get(5).position());
        assertEquals(0.0, ratios.get(5).tumorGCRatio(), 0.01);
        assertEquals(0.0, ratios.get(5).tumorReadDepth(), 0.01);
        assertEquals(0.5, ratios.get(5).tumorGcContent(), 0.01);
    }

    @Test
    public void handleDepth0WindowWholeGenome() throws Exception
    {
        sample = "depth_0_windows";
        tumorBamFile = getBam(sample);
        regionOffset = 0;

        createStandardChr1GCFile(8000);

        runCobalt(false);

        List<CobaltRatio> ratios = tumorRatioResults.get(_1);
        assertEquals(4001, ratios.get(4).position());
        assertEquals(100.0, ratios.get(4).tumorReadDepth(), 0.01);
        assertEquals(0.52, ratios.get(4).tumorGcContent(), 0.01);

        assertEquals(5001, ratios.get(5).position());
        assertEquals(0.0, ratios.get(5).tumorGCRatio(), 0.01);
        assertEquals(0.0, ratios.get(5).tumorReadDepth(), 0.01);
        assertEquals(0.5, ratios.get(5).tumorGcContent(), 0.01);
    }

    @Test
    public void gcNormalisation() throws Exception
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
        tumorBamFile = getBam(sample);
        regionOffset = 1_000;

        createStandardChr1GCFile(60_000);
        createStandardChr1PanelFile(60_000, 1.000001);

        runCobalt();

        // Normalisation calculations
        // Assign each read window to a GC bucket. Calculate the bucket median depths and smooth these median values
        // gc buckets:
        // bucket median count
        // 40 10 5
        // 41 11 5
        // 42 12 5
        // ...
        // 50 20 5
        //
        // These buckets get smoothed to:
        // 39 10/3
        // 40 21/3
        // 41 33/11
        // ...
        // 49 57/3
        // 50 39/3
        // 51 20/3

        // GC-normalise by multiplying each window's read depth by the smoothed
        // median read depth for windows belonging to the same GC bucket.
        // Note that the median depths are 0.40, 0.41 11, 0.42 12,..., 0.49 19
        // 1001, 2001, 3001, 4001, 5001: divide by 7.0
        // 6001, 7001, 8001, 9001, 10001: 11 -> divide by 11.0
        // ...

        // Finally, normalise by dividing by the mean of the non-negative values
        // so that the mean of the result is 1.0.
        // The values are:
        // (10, 10, 10, 20, 50)/7, (11, 11, 11, 21, 51)/11, (12, 12, 12, 22, 52)/12, ..., (19, 19, 19, 29, 59)/19, (20, 20, 20, 30, 60)/13
        // Total = 100/7.0 + 105/11 + 110/12 + ... + 145/19 = 101.76275
        // Count = 55, so mean = 1.850231745
        double mean = 1.850231745;

        // 1-1000 0 reads
        checkChr1TumorRatio(-1.0, 0);

        // 1001-5001 gc ratio = 1.0/7.0, depth = 10, 10, 10, 20, 50
        checkChr1TumorRatio(10.0 / 7.0 / mean, 1, 2, 3);
        checkChr1TumorRatio(20.0 / 7.0 / mean, 4);
        checkChr1TumorRatio(50.0 / 7.0 / mean, 5);

        // 6001-11_000 gc 0.41 and read depths 11, 11, 11, 21, 51
        checkChr1TumorRatio(1.0 / mean, 6, 7, 8);
        checkChr1TumorRatio(21.0 / 11.0 / mean, 9);
        checkChr1TumorRatio(51.0 / 11.0 / mean, 10);

        // 12_001-16_000 gc 0.42 and read depths 12, 12, 12, 22, 52
        checkChr1TumorRatio(1.0 / mean, 11, 12, 13);
        checkChr1TumorRatio(22.0 / 12.0 / mean, 14);
        checkChr1TumorRatio(52.0 / 12.0 / mean, 15);

        // 46_001-51_000 gc 0.49 and read depths 19, 19, 19, 29, 59
        checkChr1TumorRatio(1.0 / mean, 46, 47, 48);
        checkChr1TumorRatio(29.0 / 19.0 / mean, 49);
        checkChr1TumorRatio(59.0 / 19.0 / mean, 50);

        // 52_001-56_000 gc 0.50 and read depths 20, 20, 20, 30, 60
        checkChr1TumorRatio(20.0 / 13.0 / mean, 51, 52, 53);
        checkChr1TumorRatio(30.0 / 13.0 / mean, 54);
        checkChr1TumorRatio(60.0 / 13.0 / mean, 55);

        // 56_001-59_000 all have depth 0.
        checkChr1TumorRatio(0.0, 56, 57, 58, 59);
    }

    @Test
    public void segmentation() throws Exception
    {
        // chr1 and chr2 both of length 101_000
        // 1:1001-40_000 depth approximately 10
        // 1:40_001-70_000 depth approximately 100
        // 1:70_001-100_000 depth approximately 10
        // 2:1001-40_000 depth approximately 100
        // 2:40_001-70_000 depth approximately 10
        // 2:70_001-100_000 depth approximately 100
        // GC ratio 0.5 throughout.
        // Each window is depth +-5% (random) of the value it approximates.
        // If we have constant read depths then the segmentation algorithm assigns a cost of 0
        // to segment creation, with the result that many segments are created instead of the
        // few expected ones. (The segment cost is associated with the variance of the values being segmented.)
        sample = "multiple_segments";
        tumorBamFile = getBam(sample);
        regionOffset = 0;

        createStandardMultiChromosomeGCFile(100_000, _1, _2);
        createStandardMultipleChromosomePanelFile(100_000, 1.0001, _1, _2);
        runCobalt();

        String segmentsFile = PCFFile.generateRatioFilename(outputDir.getAbsolutePath(), sample);

        ListMultimap<Chromosome, PCFPosition> pcfData = PCFFile.readPositions(1000, PCFSource.TUMOR_BAF, segmentsFile);
        assertEquals(2, pcfData.keySet().size());
        List<PCFPosition> chr1Positions = pcfData.asMap().get(_1).stream().toList();
        // The R program that does segmentation puts a spurious 1-window segment
        // at the start of each chromosome.
        assertEquals(6, chr1Positions.size());
        assertEquals(1, chr1Positions.get(0).Position);
        assertEquals(1001, chr1Positions.get(1).Position);
        assertEquals(2001, chr1Positions.get(2).Position);
        assertEquals(41001, chr1Positions.get(3).Position);
        assertEquals(71001, chr1Positions.get(4).Position);
        assertEquals(101001, chr1Positions.get(5).Position);
    }

    @Test
    public void segmentationNew() throws Exception
    {
        sample = "multiple_segments";
        tumorBamFile = getBam(sample);
        regionOffset = 0;

        createStandardMultiChromosomeGCFile(100_000, _1, _2);
        createStandardMultipleChromosomePanelFile(100_000, 1.0001, _1, _2);
        runCobalt(false, true);

        String segmentsFileName = PCFFile.generateCobaltPcfFilename(outputDir.getAbsolutePath(), sample);
        ListMultimap<Chromosome, CobaltSegment> regions = PCFFile.readCobaltPcfFile(segmentsFileName);
        assertEquals(6, regions.size());
        List<CobaltSegment> chr1Regions = regions.get(_1);
        assertEquals(1, chr1Regions.get(0).start());
        assertEquals("1", chr1Regions.get(0).chromosome());
        assertEquals(41000, chr1Regions.get(0).end());
        assertEquals("1", chr1Regions.get(1).chromosome());
        assertEquals(41001, chr1Regions.get(1).start());
        assertEquals(71000, chr1Regions.get(1).end());
        assertEquals("1", chr1Regions.get(2).chromosome());
        assertEquals(71001, chr1Regions.get(2).start());
        assertEquals(101000, chr1Regions.get(2).end());
        assertEquals("2", regions.get(_2).get(0).chromosome());
        assertEquals(1, regions.get(_2).get(0).start());
        assertEquals(41000, regions.get(_2).get(0).end());
    }

    @Test
    public void segmentationTumorGermline() throws Exception
    {
        setupForSingleWindowBamTumorAndGermline();
        runCobalt(false, true);

        String tumorSegmentsFileName = PCFFile.generateCobaltPcfFilename(outputDir.getAbsolutePath(), sample);
        ListMultimap<Chromosome, CobaltSegment> regions = PCFFile.readCobaltPcfFile(tumorSegmentsFileName);
        assertEquals(1, regions.size());

        String referenceSegmentsFileName = PCFFile.generateCobaltPcfFilename(outputDir.getAbsolutePath(), referenceSample);
        ListMultimap<Chromosome, CobaltSegment> referenceRegions = PCFFile.readCobaltPcfFile(referenceSegmentsFileName);
        assertEquals(1, referenceRegions.size());
    }

    @Test
    public void segmentationTumorOnly() throws Exception
    {
        setupForSingleWindowBamTumorOnly();
        runCobalt(false, true);

        String tumorSegmentsFileName = PCFFile.generateCobaltPcfFilename(outputDir.getAbsolutePath(), sample);
        ListMultimap<Chromosome, CobaltSegment> regions = PCFFile.readCobaltPcfFile(tumorSegmentsFileName);
        assertEquals(1, regions.size());

        String referenceSegmentsFileName = PCFFile.generateCobaltPcfFilename(outputDir.getAbsolutePath(), referenceSample);
        assertFalse(new File(referenceSegmentsFileName).exists());
    }

    @Test
    public void segmentationGermlineOnly() throws Exception
    {
        setupForSingleWindowBamGermlineOnly();
        runCobalt(false, true);

        String tumorSegmentsFileName = PCFFile.generateCobaltPcfFilename(outputDir.getAbsolutePath(), sample);
        assertFalse(new File(tumorSegmentsFileName).exists());

        String referenceSegmentsFileName = PCFFile.generateCobaltPcfFilename(outputDir.getAbsolutePath(), referenceSample);
        ListMultimap<Chromosome, CobaltSegment> referenceRegions = PCFFile.readCobaltPcfFile(referenceSegmentsFileName);
        assertEquals(1, referenceRegions.size());
    }

    @Test
    public void wholeGenomeFemale() throws Exception
    {
        // The bam has reads for chr1, chr2 and chrX.
        // Each chromosome has length 3000.
        // 1:1-1000 depth 100
        // 1:1001-2000 depth 110
        // 1:2001-3000 depth 0
        // 2:1-1000 depth 110
        // 2:1001-2000 depth 100
        // 2:2001-3000 depth 0
        // X:1-1000 depth 120
        // X:1001-2001 depth 120
        // GC ratio 0.5 throughout.
        sample = "female";
        tumorBamFile = getBam(sample);
        regionOffset = 0;

        createStandardMultiChromosomeGCFile(3_000, _1, _2, _X);
        runCobalt(false);

        // There is only one GC bucket, and it has median read depth 105
        // (the sex chromosomes don't contribute to this median).
        // Because of GC bucket smoothing, this median gets divided by 3.
        double smoothedGcNormalisation = 105.0 / 3;
        List<CobaltRatio> ratios1 = tumorRatioResults.get(_1);
        assertEquals(3, ratios1.size());
        assertEquals(100.0 / smoothedGcNormalisation, ratios1.get(0).tumorGCRatio(), 0.01);
        assertEquals(110.0 / smoothedGcNormalisation, ratios1.get(1).tumorGCRatio(), 0.01);
        assertEquals(0.0, ratios1.get(2).tumorGCRatio(), 0.01);
        assertEquals(-1.0, ratios1.get(0).referenceGCRatio(), 0.01);
        assertEquals(-1.0, ratios1.get(0).referenceGCDiploidRatio(), 0.01);
        assertEquals(-1.0, ratios1.get(0).referenceGCRatio(), 0.01);
        assertEquals(-1.0, ratios1.get(0).referenceGCDiploidRatio(), 0.01);
        assertEquals(-1.0, ratios1.get(0).referenceGCRatio(), 0.01);
        assertEquals(-1.0, ratios1.get(0).referenceGCDiploidRatio(), 0.01);

        List<CobaltRatio> ratios2 = tumorRatioResults.get(_2);
        assertEquals(3, ratios2.size());
        assertEquals(110.0 / smoothedGcNormalisation, ratios2.get(0).tumorGCRatio(), 0.01);
        assertEquals(100.0 / smoothedGcNormalisation, ratios2.get(1).tumorGCRatio(), 0.01);
        assertEquals(0.0, ratios2.get(2).tumorGCRatio(), 0.01);

        List<CobaltRatio> ratiosX = tumorRatioResults.get(_X);
        assertEquals(3, ratiosX.size());
        assertEquals(120.0 / smoothedGcNormalisation, ratiosX.get(0).tumorGCRatio(), 0.01);
        assertEquals(120.0 / smoothedGcNormalisation, ratiosX.get(1).tumorGCRatio(), 0.01);
        assertEquals(0.0, ratiosX.get(2).tumorGCRatio(), 0.01);
    }

    @Test
    public void wholeGenomeMale() throws Exception
    {
        // The bam has reads for chr1, chr2, chrX and chrY.
        // Each chromosome has length 3000.
        // 1:1-1000 depth 100
        // 1:1001-2001 depth 110
        // 1:2001-3001 depth 104
        // 2:1-1000 depth 110
        // 2:1001-2001 depth 100
        // 2:2001-3001 depth 104
        // X:1-1000 depth 60
        // X:1001-2001 depth 60
        // X:2001-3001 depth 50
        // Y:1-1000 depth 40
        // Y:1001-2001 depth 40
        // Y:2001-3001 depth 50
        // GC ratio 0.5 throughout.
        sample = "male";
        tumorBamFile = getBam(sample);
        regionOffset = 0;

        createStandardMultiChromosomeGCFile(3_000, _1, _2, _X, _Y);
        runCobalt(false);

        // There is only one GC bucket, and it has median read depth 104
        // (the sex chromosomes don't contribute to this median).
        // Because of GC bucket smoothing, this median gets divided by 3.
        double normalisationFactor = 104.0 / 3;
        // In whole genome mode the ratios are further normalised by the
        // median/mean read depth for allosome reads (100, 110, 104, 110, 100, 104)
        normalisationFactor /= 104.0 / ((100.0 + 110.0 + 104.0) / 3.0);
        List<CobaltRatio> ratios1 = tumorRatioResults.get(_1);
        assertEquals(3, ratios1.size());
        assertEquals(100.0 / normalisationFactor, ratios1.get(0).tumorGCRatio(), 0.01);
        assertEquals(110.0 / normalisationFactor, ratios1.get(1).tumorGCRatio(), 0.01);
        assertEquals(104.0 / normalisationFactor, ratios1.get(2).tumorGCRatio(), 0.01);
        assertEquals(-1.0, ratios1.get(0).referenceGCRatio(), 0.01);
        assertEquals(-1.0, ratios1.get(0).referenceGCDiploidRatio(), 0.01);
        assertEquals(-1.0, ratios1.get(0).referenceGCRatio(), 0.01);
        assertEquals(-1.0, ratios1.get(0).referenceGCDiploidRatio(), 0.01);
        assertEquals(-1.0, ratios1.get(0).referenceGCRatio(), 0.01);
        assertEquals(-1.0, ratios1.get(0).referenceGCDiploidRatio(), 0.01);

        List<CobaltRatio> ratios2 = tumorRatioResults.get(_2);
        assertEquals(3, ratios2.size());
        assertEquals(110.0 / normalisationFactor, ratios2.get(0).tumorGCRatio(), 0.01);
        assertEquals(100.0 / normalisationFactor, ratios2.get(1).tumorGCRatio(), 0.01);
        assertEquals(104.0 / normalisationFactor, ratios2.get(2).tumorGCRatio(), 0.01);

        List<CobaltRatio> ratiosX = tumorRatioResults.get(_X);
        assertEquals(3, ratiosX.size());
        assertEquals(60.0 / normalisationFactor, ratiosX.get(0).tumorGCRatio(), 0.01);
        assertEquals(60.0 / normalisationFactor, ratiosX.get(1).tumorGCRatio(), 0.01);
        assertEquals(50.0 / normalisationFactor, ratiosX.get(2).tumorGCRatio(), 0.01);

        List<CobaltRatio> ratiosY = tumorRatioResults.get(_Y);
        assertEquals(3, ratiosY.size());
        assertEquals(40.0 / normalisationFactor, ratiosY.get(0).tumorGCRatio(), 0.01);
        assertEquals(40.0 / normalisationFactor, ratiosY.get(1).tumorGCRatio(), 0.01);
        assertEquals(50.0 / normalisationFactor, ratiosY.get(2).tumorGCRatio(), 0.01);
    }

    @Test
    public void wholeGenomeDepth0Region() throws Exception
    {
        // The tumor bam has reads for chr1, chr2 and chrX.
        // Each chromosome has length 10_000.
        // 1-3000 depth 10 for each chromosome
        // 3001-7000 depth 0 for each chromosome
        // 7001-10_000 depth 12 for each chromosome
        // The reference bam has reads for the same chromosomes and length 10_000 for each chromosome.
        // 1-2000 depth 5 for each chromosome
        // 2001-6000 depth 0 for each chromosome
        // 6001-10_000 depth 6 for each chromosome
        // GC ratio 0.5 throughout, for both tumor and reference.

        sample = "depth_0_region";
        referenceBamFile = getBam(sample + "_R");
        tumorBamFile = getBam(sample + "_T");
        regionOffset = 0;
        createStandardMultiChromosomeGCFile(10_000, _1, _2, _X);

        runCobalt(false);
        // Reference read depths are (5*2, 0*4, 6*4)*3. Median of the non-zero values is 6.0.
        // Tumor read depths are (10*3, 0*4, 12*3)*3. Mean of the non-zero values is 11.0.
        // These values get divided by three during GC bucket smoothing.
        double refGCMean = 6.0 / 3.0;
        double tumorGCMean = 11.0 / 3.0;
        // In whole genome mode the tumor and reference ratios are all divided by a normalisation
        // factor which is the median/mean of positive read depths. For this tumor this ratio is 1.0.
        // For this reference sample the value is 34.0/6.0/6.0
        double meanMedianNormalisationRef = 34.0 / 36.0;
        refGCMean = refGCMean * meanMedianNormalisationRef;
        List<CobaltRatio> ratios1 = tumorRatioResults.get(_1);
        assertEquals(10, ratios1.size());
        assertEquals(5.0, ratios1.get(0).referenceReadDepth(), 0.01);
        assertEquals(10.0, ratios1.get(0).tumorReadDepth(), 0.01);
        assertEquals(5.0 / refGCMean, ratios1.get(0).referenceGCRatio(), 0.01);
        assertEquals(5.0 / refGCMean, ratios1.get(0).referenceGCDiploidRatio(), 0.01);
        assertEquals(10.0 / tumorGCMean, ratios1.get(0).tumorGCRatio(), 0.01);

        assertEquals(5.0, ratios1.get(1).referenceReadDepth(), 0.01);
        assertEquals(10.0, ratios1.get(1).tumorReadDepth(), 0.01);
        assertEquals(5.0 / refGCMean, ratios1.get(1).referenceGCRatio(), 0.01);
        assertEquals(5.0 / refGCMean, ratios1.get(1).referenceGCDiploidRatio(), 0.01);
        assertEquals(10.0 / tumorGCMean, ratios1.get(1).tumorGCRatio(), 0.01);

        assertEquals(0.0, ratios1.get(2).referenceReadDepth(), 0.01);
        assertEquals(10.0, ratios1.get(2).tumorReadDepth(), 0.01);
        assertEquals(0.0, ratios1.get(2).referenceGCRatio(), 0.01);
        assertEquals(0.0, ratios1.get(2).referenceGCDiploidRatio(), 0.01);
        assertEquals(10.0 / tumorGCMean, ratios1.get(2).tumorGCRatio(), 0.01);

        assertEquals(0.0, ratios1.get(3).referenceReadDepth(), 0.01);
        assertEquals(0.0, ratios1.get(3).tumorReadDepth(), 0.01);
        assertEquals(0.0, ratios1.get(3).referenceGCRatio(), 0.01);
        assertEquals(0.0, ratios1.get(3).referenceGCDiploidRatio(), 0.01);
        assertEquals(0.0, ratios1.get(3).tumorGCRatio(), 0.01);

        assertEquals(0.0, ratios1.get(4).referenceReadDepth(), 0.01);
        assertEquals(0.0, ratios1.get(4).tumorReadDepth(), 0.01);
        assertEquals(0.0, ratios1.get(4).referenceGCRatio(), 0.01);
        assertEquals(0.0, ratios1.get(4).referenceGCDiploidRatio(), 0.01);
        assertEquals(0.0, ratios1.get(4).tumorGCRatio(), 0.01);

        assertEquals(0.0, ratios1.get(5).referenceReadDepth(), 0.01);
        assertEquals(0.0, ratios1.get(5).tumorReadDepth(), 0.01);
        assertEquals(0.0, ratios1.get(5).referenceGCRatio(), 0.01);
        assertEquals(0.0, ratios1.get(5).referenceGCDiploidRatio(), 0.01);
        assertEquals(0.0, ratios1.get(5).tumorGCRatio(), 0.01);

        assertEquals(6.0, ratios1.get(6).referenceReadDepth(), 0.01);
        assertEquals(0.0, ratios1.get(6).tumorReadDepth(), 0.01);
        assertEquals(6.0 / refGCMean, ratios1.get(6).referenceGCRatio(), 0.01);
        assertEquals(6.0 / refGCMean, ratios1.get(6).referenceGCDiploidRatio(), 0.01);
        assertEquals(0.0, ratios1.get(6).tumorGCRatio(), 0.01);

        assertEquals(6.0, ratios1.get(7).referenceReadDepth(), 0.01);
        assertEquals(12.0, ratios1.get(7).tumorReadDepth(), 0.01);
        assertEquals(6.0 / refGCMean, ratios1.get(7).referenceGCRatio(), 0.01);
        assertEquals(6.0 / refGCMean, ratios1.get(7).referenceGCDiploidRatio(), 0.01);
        assertEquals(12.0 / tumorGCMean, ratios1.get(7).tumorGCRatio(), 0.01);
    }

    @Test
    public void filterOutWindowsThatIntersectPseudoIgRegions() throws Exception
    {
        // The bam has reads for chr15 and chr16 only.
        // The bam registers the v38 lengths of these chromosomes, 101991189 and 90338345 respectively.
        // chr15: reads in 19_900_000-22_200_000, depth 10
        // chr16: 31_900_000-34_100_000, depth 10
        /*
        The pseudo-gene regions on chr15 and the corresponding masked-out ratios are:
        15	19964674	19964794	19964001-19965000
        15	19972790	19972910	19972001-19973000
        15	19987664	19987784	19987001-19988000
        15	22160439	22160559	22160001-22161000
        15	22178113	22178233	22178001-22179000
        15	22184975	22185095	22184001-22185000, 22185001-22186000
        15	22194891	22195011	22194001-22195000, 22195001-22196000

        The chr16 regions are:
        16	31962255	31962375	31962001-31963000
        16	32052154	32052274	32052001-32053000
        16	32066232	32066352	32066001-32067000
        16	32848026	32848146	32848001-32849000
        16	32903768	32903888	32903001-32904000
        16	32915408	32915528	32915001-32916000
        16	32995377	32995497	32995001-32996000
        16	33009496	33009616	33009001-33010000
        16	33803089	33803209	33803001-33804000
        16	33827533	33827653	33827001-33828000
        16	33844788	33844908	33844001-33845000
        16	33938345	33938465	33938001-33939000
        16	33949974	33950094	33949001-33950000, 33950001-33951000
        16	34013397	34013517	34013001-34014000
         */

        sample = "pseudo_ig_regions";
        referenceBamFile = getBam(sample);
        tumorBamFile = getBam(sample);
        regionOffset = 0;
        int chr15Length = 101_991_189;
        int chr16Length = 90_338_345;
        Map<HumanChromosome, Integer> lengthsMap = new HashMap<>();
        lengthsMap.put(_15, rounded1000(chr15Length));
        lengthsMap.put(_16, rounded1000(chr16Length));
        createStandardMultiChromosomeGCFile(lengthsMap);

        runCobalt(false);

        // We check that for each of the expected masked regions, the ratio is -1.0,
        // but for the surrounding windows, the ratio is 1.0.
        // A more ambitious test would be to check every single window.
        // This is complicated by some regions of 'N' bases between the Ig
        // pseudo gene regions on chr 15.
        List<Integer> chr15Exclusions = new ArrayList<>();
        chr15Exclusions.add(19964001);
        chr15Exclusions.add(19972001);
        chr15Exclusions.add(19987001);
        chr15Exclusions.add(22160001);
        chr15Exclusions.add(22178001);
        chr15Exclusions.add(22184001);
        chr15Exclusions.add(22194001);
        chr15Exclusions.add(22185001);
        chr15Exclusions.add(22195001);
        checkWindowsAreMaskedButNeighoursAreNot(chr15Exclusions, _15);

        List<Integer> chr16Exclusions = new ArrayList<>();
        chr16Exclusions.add(31962001);
        chr16Exclusions.add(32052001);
        chr16Exclusions.add(32066001);
        chr16Exclusions.add(32848001);
        chr16Exclusions.add(32903001);
        chr16Exclusions.add(32915001);
        chr16Exclusions.add(32995001);
        chr16Exclusions.add(33009001);
        chr16Exclusions.add(33803001);
        chr16Exclusions.add(33827001);
        chr16Exclusions.add(33844001);
        chr16Exclusions.add(33938001);
        chr16Exclusions.add(33949001);
        chr16Exclusions.add(34013001);
        chr16Exclusions.add(33950001);
        checkWindowsAreMaskedButNeighoursAreNot(chr16Exclusions, _16);
    }

    @Test
    public void regionsOfExtremeGCAreFilteredOutInWholeGenomeMode() throws Exception
    {
        // 1 chr of length 101_000
        // 1:1-1000 depth 10, gc ratio = 0.00
        // 1:1001-2000 depth 10, gc ratio = 0.01
        // 1:2001-3000 depth 10, gc ratio = 0.02
        // etc
        sample = "increasing_gc_per_window_depth_10";
        tumorBamFile = getBam(sample);
        regionOffset = 1000;

        createStandardChr1GCFile(101_000);

        runCobalt(false);
        List<CobaltRatio> ratios = tumorRatioResults.get(_1);
        for(int i = 0; i < 25; i++)
        {
            assertEquals(-1.0, ratios.get(i).tumorGCRatio(), 0.01);
        }
        for(int i = 25; i < 68; i++)
        {
            assertEquals(1.0, ratios.get(i).tumorGCRatio(), 0.01);
        }
        for(int i = 68; i < 100; i++)
        {
            assertEquals(-1.0, ratios.get(i).tumorGCRatio(), 0.01);
        }
    }

    private void checkWindowsAreMaskedButNeighoursAreNot(final List<Integer> exclusions, final HumanChromosome humanChromosome)
    {
        for(Integer position : exclusions)
        {
            // The window at this position should be masked.
            int windowStart = rounded1000(position) + 1;
            assertEquals(-1.0, retrieveRatio(humanChromosome, windowStart).tumorGCRatio(), 0.01);

            // The previous window should not be masked, unless it in the masked set.
            int previousWindowStart = windowStart - 1000;
            double expectedRatio = exclusions.contains(previousWindowStart) ? -1.0 : 1.0;
            assertEquals(expectedRatio, retrieveRatio(humanChromosome, previousWindowStart).tumorGCRatio(), 0.01);

            // The next window should not be masked, unless it is in the masked set.
            int nextWindowStart = windowStart - 1000;
            expectedRatio = exclusions.contains(nextWindowStart) ? -1.0 : 1.0;
            assertEquals(expectedRatio, retrieveRatio(humanChromosome, nextWindowStart).tumorGCRatio(), 0.01);
        }
    }

    private static int rounded1000(final Integer position)
    {
        return (position / 1000) * 1000;
    }

    private CobaltRatio retrieveRatio(HumanChromosome chromosome, int position)
    {
        return tumorRatioResults.get(chromosome)
                .stream()
                .filter(t -> t.position() == position)
                .findFirst()
                .orElse(null);
    }

    private void checkChr1TumorRatio(double expected, int... indices)
    {
        List<CobaltRatio> ratios = tumorRatioResults.get(_1);
        for(final int index : indices)
        {
            assertEquals(expected, ratios.get(index).tumorGCRatio(), 0.01);
        }
    }

    private void setupForSingleWindowBamTumorOnly() throws IOException
    {
        setupForSingleWindowBamTumorAndGermline();
        referenceBamFile = null;
        referenceSample = null;
    }

    private void setupForSingleWindowBamGermlineOnly() throws IOException
    {
        setupForSingleWindowBamTumorAndGermline();
        tumorBamFile = null;
    }

    private void setupForSingleWindowBamTumorAndGermline() throws IOException
    {
        sample = "one_window";
        referenceSample = sample + "_reference";
        referenceBamFile = getBam(sample);
        tumorBamFile = getBam(sample);
        regionOffset = 0;

        createStandardChr1GCFile(3_000);
        createStandardChr1PanelFile(3_000, 1.000001);
    }

    private void setupForThreeWindowBamTumorOnly() throws IOException
    {
        // 1 chr of length 5000
        // 1:1001-4000 depth 100
        sample = "three_windows";
        tumorBamFile = getBam(sample);
        regionOffset = 0;

        createStandardChr1GCFile(4_000);
        createStandardChr1PanelFile(4_000, 1.0000001);
    }

    private void createStandardMultipleChromosomePanelFile(final int length, final double relativeEnrichment,
            HumanChromosome... chromosomes) throws IOException
    {
        panelNormalisation = new File(tempDir, "ThePanel.tsv");
        PanelFileWriter panelWriter = new PanelFileWriter();
        for(HumanChromosome chr : chromosomes)
        {
            panelWriter.addSection(new PanelFileSection(chr, regionOffset, regionOffset + length, relativeEnrichment));
        }
        panelWriter.write(panelNormalisation);
    }

    private void createStandardChr1PanelFile(final int length, final double relativeEnrichment) throws IOException
    {
        createStandardMultipleChromosomePanelFile(length, relativeEnrichment, _1);
    }

    private void createStandardMultiChromosomeGCFile(final int length, HumanChromosome... chromosomes) throws IOException
    {
        Map<HumanChromosome, Integer> lengthsMap = new HashMap<>();
        for(HumanChromosome chr : chromosomes)
        {
            lengthsMap.put(chr, length);
        }
        createStandardMultiChromosomeGCFile(lengthsMap);
    }

    private void createStandardMultiChromosomeGCFile(Map<HumanChromosome, Integer> lengthsMap) throws IOException
    {
        gcProfile = new File(tempDir, "GC_profile.1000bp.38.cnp");
        GcProfilesUtilities gcFileWriter = new GcProfilesUtilities();
        for(HumanChromosome chr : lengthsMap.keySet())
        {
            int length = lengthsMap.get(chr);
            gcFileWriter.addSection(new ConstantGcFileSection(chr, regionOffset, regionOffset + length, 0.5));
        }
        gcFileWriter.write(gcProfile);
    }

    private void createStandardChr1GCFile(final int length) throws IOException
    {
        createStandardMultiChromosomeGCFile(length, _1);
    }

    private File getBam(String sample)
    {
        return Path.of("src", "test", "resources", "bam", sample + ".sorted.bam").toFile();
    }

    private void runCobalt() throws Exception
    {
        runCobalt(true);
    }

    private void runCobalt(boolean targeted) throws Exception
    {
        runCobalt(targeted, false);
    }

    private void runCobalt(boolean targeted, boolean useNewSegmenter) throws Exception
    {
        int argCount = 12;
        if(targeted)
        {
            argCount += 2;
        }
        if(useNewSegmenter)
        {
            argCount += 1;
        }
        if(tumorBamFile != null && referenceBamFile != null)
        {
            argCount += 4;
        }
        if(diploidBedFile != null)
        {
            argCount += 2;
        }
        String[] args = new String[argCount];
        int index = 0;
        args[index++] = String.format("-%s", GC_PROFILE);
        args[index++] = String.format("%s", gcProfile.getAbsolutePath());
        args[index++] = String.format("-%s", PCF_GAMMA);
        args[index++] = String.format("%d", 50);
        args[index++] = String.format("-%s", OUTPUT_DIR);
        args[index++] = String.format("%s", outputDir.getAbsolutePath());
        args[index++] = String.format("-%s", REF_GENOME_VERSION);
        args[index++] = String.format("%s", "38");
        if(tumorBamFile != null)
        {
            args[index++] = String.format("-%s", TUMOR);
            args[index++] = String.format("%s", sample);
            args[index++] = String.format("-%s", TUMOR_BAM);
            args[index++] = String.format("%s", tumorBamFile.getAbsolutePath());
        }
        if(referenceBamFile != null)
        {
            args[index++] = String.format("-%s", REFERENCE);
            args[index++] = String.format("%s", referenceSample);
            args[index++] = String.format("-%s", REFERENCE_BAM);
            args[index++] = String.format("%s", referenceBamFile.getAbsolutePath());
        }
        if(targeted)
        {
            args[index++] = String.format("-%s", TARGET_REGION_NORM_FILE);
            args[index++] = String.format("%s", panelNormalisation.getAbsolutePath());
        }
        if(useNewSegmenter)
        {
            args[index++] = String.format("-%s", USE_NEW_SEGMENTER);
        }
        if(diploidBedFile != null)
        {
            args[index++] = String.format("-%s", TUMOR_ONLY_DIPLOID_BED);
            args[index] = String.format("%s", diploidBedFile.getAbsolutePath());
        }

        CobaltApplication.main(args);

        File ratioFile;
        if(tumorBamFile != null)
        {
            ratioFile = new File(outputDir, sample + ".cobalt.ratio.tsv.gz");
        }
        else
        {
            ratioFile = new File(outputDir, referenceSample + ".cobalt.ratio.tsv.gz");
        }
        assertTrue(ratioFile.exists());
        assertTrue(ratioFile.isFile());
        RawCobaltRatioFile rawResultsFile = new RawCobaltRatioFile(ratioFile.getAbsolutePath());
        List<RawCobaltRatio> rawRatios = rawResultsFile.read();
        tumorRatioResults = new HashMap<>();
        for(RawCobaltRatio rawRatio : rawRatios)
        {
            Chromosome chromosome = HumanChromosome.fromString(rawRatio.chromosome());
            List<CobaltRatio> chrResults;
            if(tumorRatioResults.containsKey(chromosome))
            {
                chrResults = tumorRatioResults.get(chromosome);
            }
            else
            {
                chrResults = new ArrayList<>();
                tumorRatioResults.put(chromosome, chrResults);
            }
            chrResults.add(rawRatio.toCobaltRatio());
        }

        File medianRatiosFile = new File(CobaltMedianRatioFile.generateFilename(outputDir.getAbsolutePath(), sample));
        if(medianRatiosFile.exists())
        {
            medianRatioResults = CobaltMedianRatioFile.read(medianRatiosFile.getAbsolutePath());
        }
        File referenceMedianRatiosFile = new File(CobaltMedianRatioFile.generateFilename(outputDir.getAbsolutePath(), referenceSample));
        if(referenceMedianRatiosFile.exists())
        {
            referenceMedianRatioResults = CobaltMedianRatioFile.read(referenceMedianRatiosFile.getAbsolutePath());
        }
    }
}