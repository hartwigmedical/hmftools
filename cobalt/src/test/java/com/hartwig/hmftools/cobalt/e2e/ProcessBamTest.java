package com.hartwig.hmftools.cobalt.e2e;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.CobaltApplication;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.purple.Gender;

import org.apache.commons.io.FileUtils;
import org.junit.Assert;
import org.junit.Test;

public class ProcessBamTest
{
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
        args[index++] = String.format("-%s", "target_region_norm_file");
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
