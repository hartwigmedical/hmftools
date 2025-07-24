package com.hartwig.hmftools.cobalt.e2e;

import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;

import java.io.File;

import com.hartwig.hmftools.cobalt.CobaltApplication;

import org.apache.commons.io.FileUtils;
import org.junit.Test;

public class NoDeletionsTest
{
    @Test
    public void noDeletions() throws Exception
    {
        File tempDir = new File("/Users/timlavers/work/junk/rubbish");
        File bamFile = new File(tempDir, "T2.sorted.bam");
        int regionOffset = 10_000_000;

        File gcProfile = new File(tempDir, "GC_profile.1000bp.38.cnp");
        new GcProfilesUtilities().writeDefaultGcProfile(regionOffset + 10_000, regionOffset + 40_000, gcProfile);

        File panelNormalisation = new File(tempDir, "ThePanel.tsv");
        new PanelFileUtilities().writePanelNormalisationFile(regionOffset + 10_000, regionOffset + 40_000, panelNormalisation);

        File refGenome = new File("/Users/timlavers/work/data/reference_genome_no_alts/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna");

        File outputDir = new File(tempDir, "output");
        outputDir.mkdirs();
        FileUtils.cleanDirectory(outputDir);

        File debugDir = new File(tempDir, "debug");
        debugDir.mkdirs();
        FileUtils.cleanDirectory(debugDir);

        String[] args = new String[14];
        int index = 0;
        args[index++] = String.format("-%s", TUMOR);
        args[index++] = String.format("%s", "T1");
        args[index++] = String.format("-%s", TUMOR_BAM);
        args[index++] = String.format("%s", bamFile.getAbsolutePath());
        args[index++] = String.format("-%s", GC_PROFILE);
        args[index++] = String.format("%s", gcProfile.getAbsolutePath());
        args[index++] = String.format("-%s", "target_region_norm_file");
        args[index++] = String.format("%s", panelNormalisation.getAbsolutePath());
        args[index++] = String.format("-%s", REF_GENOME);
        args[index++] = String.format("%s", refGenome.getAbsolutePath());
        args[index++] = String.format("-%s", "pcf_gamma");
        args[index++] = String.format("%d", 1);
        args[index++] = String.format("-%s", OUTPUT_DIR);
        args[index] = String.format("%s", outputDir.getAbsolutePath());

        CobaltApplication.main(args);
    }
}
