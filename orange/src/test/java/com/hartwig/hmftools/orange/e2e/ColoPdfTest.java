package com.hartwig.hmftools.orange.e2e;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CUPPA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LILAC_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_PLOT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PEACH_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_PLOT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.QSEE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_PLOT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SIGS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.orange.OrangeConfig.EXPERIMENT_TYPE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;

import com.hartwig.hmftools.orange.OrangeApplication;
import com.hartwig.hmftools.common.test.Unzipper;

import org.apache.commons.io.FileUtils;
import org.apache.pdfbox.Loader;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.PDPage;
import org.apache.pdfbox.text.PDFTextStripper;
import org.junit.Before;
import org.junit.Test;

public class ColoPdfTest
{
    String tumor = "COLO829v004T";
    String reference = "COLO829v004R";

    private File tempDir;
    private File outputDir;
    private File inputsDir;

    @Before
    public void setUp() throws IOException
    {
        //        tempDir = Files.createTempDirectory("pbt").toFile();
        tempDir = new File("/Users/timlavers/work/batches/2026/6/4/1/tests");
        outputDir = new File(tempDir, "output");
        inputsDir = new File(tempDir, "inputs");
        //noinspection ResultOfMethodCallIgnored
        outputDir.mkdirs();
        FileUtils.cleanDirectory(outputDir);
        //noinspection ResultOfMethodCallIgnored
        inputsDir.mkdirs();
        FileUtils.cleanDirectory(inputsDir);
    }

    @Test
    public void testColoPdf() throws Exception
    {
        final String dirName = "colo_swapped_images";
        File inputsZip = getInputsZip(dirName);
        Unzipper.unzipInto(inputsZip, inputsDir);
        runOrange(dirName);

        File outputFile = new File(outputDir, tumor + ".orange.pdf");
        assertTrue(outputFile.exists());

        try(PDDocument document = Loader.loadPDF(outputFile))
        {
            assertEquals(8, document.getNumberOfPages());

            // Check specific text on Page 1
            PDFTextStripper stripper = new PDFTextStripper();
            stripper.setStartPage(1);
            stripper.setEndPage(1);
            String page1Text = stripper.getText(document);
            String[] lines = page1Text.split("\n");

            assertEquals("ORANGE Report", lines[0]);
            assertEquals("SAMPLE", lines[1]);
            assertEquals(tumor, lines[2]);
            assertEquals("1/8", lines[3]);
            // First table
            assertEquals("PRIMARY TUMOR PURITY PLOIDY FIT METHOD QC", lines[4]);
            assertEquals("NOT SPECIFIED 99% (97%-100%) 3.00 (2.96-3.05) NORMAL PASS", lines[5]);

            // Second table
            assertEquals("PIPELINE VERSION GENOME VERSION SEQUENCING TYPE PIPELINE SAMPLES DATE ANALYSED", lines[6]);
            assertTrue(lines[7].contains("V38 ILLUMINA WHOLE GENOME TUMOR / NORMAL"));

            // Side-by-side tables: Driver Summary (left) and Genome Wide Biomarkers (right)
            PDPage page1 = document.getPage(0);
            PageTextRipper page1Ripper = new PageTextRipper(page1);

            // The driver summary is on the left side of the page, from about 23% to about 50% of the way down.
            String[] driverSummary = page1Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.23), new PositionInPage(0.5, 0.5));
            assertEquals("Driver Summary", driverSummary[0]);
            assertEquals("Somatic variant: 7 (BRAF, CDKN2A, HDAC2, TERT)", driverSummary[1]);
            assertEquals("Somatic copy number: 1 (PTEN)", driverSummary[2]);
            assertEquals("Somatic disruption: 2", driverSummary[3]);
            assertEquals("Germline variant: None", driverSummary[4]);
            assertEquals("Germline copy number: None", driverSummary[5]);
            assertEquals("Germline disruption: None", driverSummary[6]);
            assertEquals("Fusion drivers: None", driverSummary[7]);
            assertEquals("Viral presence: None", driverSummary[8]);
            assertEquals("Whole genome Yes", driverSummary[9]);
            assertEquals("duplicated:", driverSummary[10]);
            assertEquals("DPYD status: *1 HOM (Normal Function)", driverSummary[11]);

            // The Genome Wide Biomarkers table is to the right of the driver summary table
            String[] genomeWideBiomarkers = page1Ripper.getLinesInRectangle(new PositionInPage(0.5, 0.23), new PositionInPage(1.0, 0.5));
            assertEquals("Genome Wide Biomarkers", genomeWideBiomarkers[0]);
            assertEquals("Microsatellite indels per Mb: 0.1 (Stable)", genomeWideBiomarkers[1]);
            assertEquals("Tumor mutations per Mb: 14.4 (High)", genomeWideBiomarkers[2]);
            assertEquals("Tumor mutational load: 198 (High)", genomeWideBiomarkers[3]);
            assertEquals("HR deficiency score: 0.0 (Proficient)", genomeWideBiomarkers[4]);
            assertEquals("LOH proportion: 15%", genomeWideBiomarkers[5]);
            assertEquals("Number of SVs: 125", genomeWideBiomarkers[6]);
            assertEquals("CUPPA cancer type: Skin: Melanoma (100%)", genomeWideBiomarkers[7]);
        }
    }

    /*
All 14 images have been rewritten. Here's the colour reference:
| File                    | Colour        | RGB             |
| ----------------------- | ------------- | --------------- |
| `copynumber.png`        | Red-orange    | (255, 69, 0)    |
| `somatic.rainfall.png`  | Gold          | (255, 215, 0)   |
| `cluster-62.003.png`    | Lime green    | (50, 205, 50)   |
| `cluster-12.001.png`    | Deep sky blue | (0, 191, 255)   |
| `somatic.png`           | Violet        | (148, 0, 211)   |
| `input.png`             | Hot pink      | (255, 105, 180) |
| `map.png`               | Dark orange   | (255, 140, 0)   |
| `purity.range.png`      | Teal          | (0, 128, 128)   |
| `segment.png`           | Crimson       | (220, 20, 60)   |
| `circos.png`            | Indigo        | (75, 0, 130)    |
| `cluster-59.007.png`    | Chartreuse    | (127, 255, 0)   |
| `somatic.clonality.png` | Dodger blue   | (30, 144, 255)  |
| `cuppa.vis.png`         | Deep pink     | (255, 20, 147)  |
| `qsee.vis.report.png`   | Saddle brown  | (139, 69, 19)   |
     */

    private void runOrange(String dirName) throws Exception
    {
        int argCount = 37;
        String[] args = new String[argCount];
        int index = 0;
        final File inputsDirectory = new File(inputsDir, dirName);
        String inputsDirStr = inputsDirectory.getAbsolutePath();
        String plotsDirStr = new File(inputsDirectory, "plot").getAbsolutePath();
        args[index++] = String.format("-%s", TUMOR);
        args[index++] = String.format("%s", tumor);
        args[index++] = String.format("-%s", REFERENCE);
        args[index++] = String.format("%s", reference);
        args[index++] = String.format("-%s", REF_GENOME_VERSION);
        args[index++] = String.format("%s", "38");
        args[index++] = String.format("-%s", LOG_DEBUG);
        args[index++] = String.format("-%s", EXPERIMENT_TYPE);
        args[index++] = String.format("%s", "WGS");

        args[index++] = String.format("-%s", LINX_GERMLINE_DIR_CFG);
        args[index++] = String.format("%s", inputsDirStr);
        args[index++] = String.format("-%s", PURPLE_DIR_CFG);
        args[index++] = String.format("%s", inputsDirStr);
        args[index++] = String.format("-%s", PURPLE_PLOT_DIR_CFG);
        args[index++] = String.format("%s", plotsDirStr);

        args[index++] = String.format("-%s", LINX_DIR_CFG);
        args[index++] = String.format("%s", inputsDirStr);
        args[index++] = String.format("-%s", LINX_PLOT_DIR_CFG);
        args[index++] = String.format("%s", plotsDirStr);

        args[index++] = String.format("-%s", QSEE_DIR_CFG);
        args[index++] = String.format("%s", inputsDirStr);

        args[index++] = String.format("-%s", LILAC_DIR_CFG);
        args[index++] = String.format("%s", inputsDirStr);

        args[index++] = String.format("-%s", SAGE_PLOT_DIR_CFG);
        args[index++] = String.format("%s", plotsDirStr);

        args[index++] = String.format("-%s", CHORD_DIR_CFG);
        args[index++] = String.format("%s", inputsDirStr);

        args[index++] = String.format("-%s", CUPPA_DIR_CFG);
        args[index++] = String.format("%s", inputsDirStr);

        args[index++] = String.format("-%s", PEACH_DIR_CFG);
        args[index++] = String.format("%s", inputsDirStr);

        args[index++] = String.format("-%s", SIGS_DIR_CFG);
        args[index++] = String.format("%s", inputsDirStr);

        args[index++] = String.format("-%s", VIRUS_DIR_CFG);
        args[index++] = String.format("%s", inputsDirStr);

        args[index++] = String.format("-%s", OUTPUT_DIR);
        args[index] = String.format("%s", outputDir.getAbsolutePath());
        OrangeApplication.main(args);
    }

    private File getInputsZip(String name)
    {
        return Path.of("src", "test", "resources", "e2e", name + ".zip").toFile();
    }
}
