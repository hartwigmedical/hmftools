package com.hartwig.hmftools.orange.e2e;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.purple.PurpleQCStatus.FAIL_CONTAMINATION;
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

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import com.hartwig.hmftools.orange.OrangeApplication;
import com.hartwig.hmftools.common.test.Unzipper;

import org.apache.commons.io.FileUtils;
import org.apache.pdfbox.Loader;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.PDPage;
import org.junit.Assert;
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
        tempDir = new File("/Users/timlavers/work/batches/2026/6/29/1");
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
    public void testColoPdfWithQcFail() throws Exception
    {
        // This test unzips the required data for generating the Orange report
        // for Colo829, overwrites the QC file, runs Orange, and checks that the
        // report has the appropriate warning message and is mostly empty.
        final String dirName = "colo_swapped_images";
        File inputsZip = getInputsZip(dirName);
        Unzipper.unzipInto(inputsZip, inputsDir);

        File filesDir = new File(inputsDir, dirName);
        File purpleQcFile = new File(filesDir, tumor + ".purple.qc");
        List<String> lines = Files.readAllLines(purpleQcFile.toPath());
        lines.set(0, "QCStatus\t" + FAIL_CONTAMINATION.name());
        Files.write(purpleQcFile.toPath(), lines);

        runOrange(dirName);

        File outputFile = new File(outputDir, tumor + ".orange.pdf");
        assertTrue(outputFile.exists());
        Assert.assertTrue(outputFile.length() > 0);

        try(PDDocument document = Loader.loadPDF(outputFile))
        {
            assertEquals(7, document.getNumberOfPages());

            PDPage page1 = document.getPage(0);
            PageRipper page1Ripper = new PageRipper(page1);
            String[] sampleTable = page1Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.08), new PositionInPage(1.0, 0.2));
            assertEquals(3, sampleTable.length);
            assertEquals("PRIMARY TUMOR PURITY PLOIDY FIT METHOD QC", sampleTable[0]);
            assertEquals("NOT SPECIFIED 99% (97%-100%) 3.00 (2.96-3.05) NORMAL FAIL_CONTAMINATION", sampleTable[1]);
            assertEquals("The QC status of this sample is fail (contamination): all presented data in this report should be interpreted with caution", sampleTable[2]);

            for(int i = 2; i <= 7; i++)
            {
                checkPageIsNA(document, i);
            }
        }
    }

    @Test
    public void testColoPdf() throws Exception
    {
        // This test unzips the required data for generating the Orange report
        // for Colo829, runs Orange, and checks the text and images of the generated PDF.
        // The original image files in the input data have been overwritten with image
        // data for images of the same size but with each image consisting of a single colour.
        // We can use the expected images sizes and colours to check that these are being
        // rendered correctly.
        final String dirName = "colo_swapped_images";
        File inputsZip = getInputsZip(dirName);
        Unzipper.unzipInto(inputsZip, inputsDir);
        runOrange(dirName);

        File outputFile = new File(outputDir, tumor + ".orange.pdf");
        assertTrue(outputFile.exists());
        Assert.assertTrue(outputFile.length() > 0);

        try(PDDocument document = Loader.loadPDF(outputFile))
        {
            //            assertEquals(9, document.getNumberOfPages());

            checkPage1(document);
            checkPage2(document);
            //            checkPage3(document);
            //            checkPage4(document);
            //            checkPage5(document);
            //            checkPage6(document);
            //            checkPage7(document);
            //            checkPage8(document);
            //            checkPage9(document);
        }
    }

    private void checkPageIsNA(PDDocument document, int page) throws IOException
    {
        PDPage pdPage = document.getPage(page - 1);
        PageRipper pageRipper = new PageRipper(pdPage);

        String[] body = pageRipper.getLinesInRectangle(new PositionInPage(0.0, 0.1), new PositionInPage(1.0, 0.4));
        assertEquals(2, body.length);
        assertEquals("NA", body[1]);
    }

    private void checkPage1(final PDDocument document) throws IOException
    {
        PDPage page1 = document.getPage(0);
        PageRipper page1Ripper = new PageRipper(page1);
        checkHeaderAndFooter(page1Ripper, 1, true);

        String[] sampleTable = page1Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.08), new PositionInPage(1.0, 0.14));
        assertEquals(2, sampleTable.length);
        assertEquals("PRIMARY TUMOR PURITY PLOIDY FIT METHOD QC", sampleTable[0]);
        assertEquals("NOT SPECIFIED 99% (97%-100%) 3.00 (2.96-3.05) NORMAL PASS", sampleTable[1]);

        String[] pieplineTable = page1Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.14), new PositionInPage(1.0, 0.19));
        assertEquals(2, pieplineTable.length);
        assertEquals("PIPELINE VERSION GENOME VERSION SEQUENCING TYPE PIPELINE SAMPLES DATE ANALYSED", pieplineTable[0]);
        assertTrue(pieplineTable[1].contains("V38 ILLUMINA WHOLE GENOME TUMOR / NORMAL"));

        // The driver summary is on the left side of the page, from about 20% to about 42% of the way down.
        String[] driverSummary = page1Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.20), new PositionInPage(0.5, 0.42));
        assertEquals("Driver Summary", driverSummary[0]);
        assertEquals("Somatic variant: 7 (BRAF, CDKN2A, HDAC2, TERT)", driverSummary[1]);
        assertEquals("Somatic copy number: 1 (PTEN)", driverSummary[2]);
        assertEquals("Somatic disruption: 2", driverSummary[3]);
        assertEquals("Germline variant: None", driverSummary[4]);
        assertEquals("Germline copy number: None", driverSummary[5]);
        assertEquals("Germline disruption: None", driverSummary[6]);
        assertEquals("Fusion drivers: None", driverSummary[7]);
        assertEquals("Viral presence: None", driverSummary[8]);
        assertEquals("Whole genome duplicated: Yes", driverSummary[9]);
        assertEquals("DPYD status: *1 HOM (Normal Function)", driverSummary[10]);

        // The Genome Wide Biomarkers table is to the right of the driver summary table
        String[] genomeWideBiomarkers = page1Ripper.getLinesInRectangle(new PositionInPage(0.5, 0.20), new PositionInPage(1.0, 0.42));
        assertEquals("Genome Wide Biomarkers", genomeWideBiomarkers[0]);
        assertEquals("Microsatellite indels per Mb: 0.1 (Stable)", genomeWideBiomarkers[1]);
        assertEquals("Tumor mutations per Mb: 14.4 (High)", genomeWideBiomarkers[2]);
        assertEquals("Tumor mutational load: 198 (High)", genomeWideBiomarkers[3]);
        assertEquals("HR deficiency score: 0.0 (Proficient)", genomeWideBiomarkers[4]);
        assertEquals("LOH proportion: 15%", genomeWideBiomarkers[5]);
        assertEquals("Number of SVs: 125", genomeWideBiomarkers[6]);
        assertEquals("CUPPA cancer type: Skin: Melanoma (100%)", genomeWideBiomarkers[7]);

        // Check the images
        assertEquals(2, page1Ripper.numberOfImages());
        final RectangleInPage topLeftCorner = new RectangleInPage(0.0, 0.0, 0.2, 0.09);
        page1Ripper.checkHasImageWithinBoundsOfGivenSize(topLeftCorner, 1234, 1200);
        final RectangleInPage lowerMiddle = new RectangleInPage(0.15, 0.4, 0.85, 0.92);
        final Color indigo = new Color(75, 0, 130);
        page1Ripper.checkHasImageWithinBoundsOfGivenSizeAndColor(lowerMiddle, 3000, 3000, indigo);
    }

    private void checkPage2(final PDDocument document) throws IOException
    {
        PDPage page2 = document.getPage(1);
        PageRipper page2Ripper = new PageRipper(page2);
        checkHeaderAndFooter(page2Ripper, 2, true);

        String[] somaticFindingsHeading = page2Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.1), new PositionInPage(1.0, 0.15));
        assertEquals(1, somaticFindingsHeading.length);
        assertEquals("Somatic Findings", somaticFindingsHeading[0]);

        String[] smallVariantsTable = page2Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.15), new PositionInPage(1.0, 0.35));
        assertEquals(10, smallVariantsTable.length);
        assertEquals("Small Variants (8)", smallVariantsTable[0]);
        assertEquals("GENE POSITION HGVS AF DEPTH COPIES HOTSPOT BIALL CLONAL DRIVER", smallVariantsTable[1]);
        assertEquals("BRAF chr7:140753336 c.1799T>A [p.Val600Glu] 0.62 237 3.8 of 6.1 Yes 2% 100% 100%", smallVariantsTable[2]);
        assertEquals("CDKN2A chr9:21971154 c.203_204delCG [p.Ala68fs] 1.00 96 2.0 of 2.0 Near 100% 100% 100%", smallVariantsTable[3]);
        assertEquals("CDKN2A chr9:21971154 c.246_247delCG [p.Gly83fs] 1.00 96 2.0 of 2.0 Near 100% 100% 100%", smallVariantsTable[4]);
        assertEquals("TERT chr5:1295113 c.-125_-124delCCinsTT 0.86 73 1.7 of 2.0 Yes 93% 100% 100%", smallVariantsTable[5]);
        assertEquals("HDAC2 chr6:113943504 c.1225C>T [p.Arg409*] 0.45 94 1.0 of 2.2 No 2% 100% 33%", smallVariantsTable[6]);
        assertEquals("GRIN2A chr16:10180333 c.79G>T [p.Ala27Ser] 0.07 121 0.2 of 3.1 No 2% 0% 14%", smallVariantsTable[7]);
        assertEquals("SF3B1 chr2:197402055 c.2153C>T [p.Pro718Leu] 0.65 115 2.0 of 3.0 No 2% 100% 14%", smallVariantsTable[8]);
        assertEquals("TP63 chr3:189886541 c.1497G>T [p.Met499Ile] 0.55 150 2.2 of 4.0 No 2% 100% 0%", smallVariantsTable[9]);

        String[] ampsAndDelsTable = page2Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.35), new PositionInPage(1.0, 0.45));
        assertEquals(3, ampsAndDelsTable.length);
        assertEquals("Amplifications and Deletions (1)", ampsAndDelsTable[0]);
        assertEquals("LOCATION GENE TYPE RANGE MIN CN MAX CN REL CN ARM CN DRIVER", ampsAndDelsTable[1]);
        assertEquals("chr10q23.31 PTEN DEL PARTIAL 0.0 2.0 0.0 2.0 HIGH", ampsAndDelsTable[2]);

        String[] fusionsTable = page2Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.45), new PositionInPage(1.0, 0.5));
        assertEquals(2, fusionsTable.length);
        assertEquals("Fusions (0)", fusionsTable[0]);
        assertEquals("NONE", fusionsTable[1]);

        String[] disruptionsTable = page2Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.5), new PositionInPage(1.0, 0.6));
        assertEquals(3, disruptionsTable.length);
        assertEquals("Disruptions (1)", disruptionsTable[0]);
        assertEquals("GENE POSITION ZYGOSITY CONTEXT TYPE JCN UNDISRUPTED CN DRIVER", disruptionsTable[1]);
        assertEquals("PTEN chr10:87940542 - chr10:87952584 HOM Intron 5 - Intron 6 DEL 2.0 0.0 LOW", disruptionsTable[2]);

        String[] virusesTable = page2Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.6), new PositionInPage(1.0, 0.65));
        assertEquals(2, virusesTable.length);
        assertEquals("Viruses (0)", virusesTable[0]);
        assertEquals("NONE", virusesTable[1]);

        String[] chrArmsTable = page2Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.65), new PositionInPage(1.0, 0.95));
        //        assertEquals(13, chrArmsTable.length);
        //        assertEquals("Arm Copy Number Aberrations", chrArmsTable[0]);
        assertEquals("CHROMOSOME ARM TYPE CN REL CN DRIVER", chrArmsTable[1]);
        assertEquals("chr1 Q GAIN 3.9 1.4 HIGH", chrArmsTable[2]);
        //        assertEquals("THE TABLE CONTINUES ON THE NEXT PAGE", chrArmsTable[12]);

        // Check the image
        assertEquals(1, page2Ripper.numberOfImages());
        final RectangleInPage topLeftCorner = new RectangleInPage(0.0, 0.0, 0.2, 0.09);
        page2Ripper.checkHasImageWithinBoundsOfGivenSize(topLeftCorner, 1234, 1200);
    }

    private void checkPage3(final PDDocument document) throws IOException
    {
        PageRipper page3Ripper = new PageRipper(document.getPage(2));
        checkHeaderAndFooter(page3Ripper, 3, true);

        String[] chrArmsTable = page3Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.1), new PositionInPage(1.0, 0.4));
        assertEquals(15, chrArmsTable.length);
        //        assertEquals("Arm Copy Number Aberrations", chrArmsTable[0]);
        assertEquals("CONTINUED FROM THE PREVIOUS PAGE", chrArmsTable[1]);
        assertEquals("CHROMOSOME ARM TYPE CN REL CN DRIVER", chrArmsTable[2]);
        assertEquals("chr7 Q GAIN 4.0 1.4 LOW", chrArmsTable[3]);
        assertEquals("chrX Q GAIN 1.9 1.3 LOW", chrArmsTable[14]);

        String[] signAllocTable = page3Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.4), new PositionInPage(1.0, 0.65));
        assertEquals(12, signAllocTable.length);
        assertEquals("Signature Allocations (12)", signAllocTable[0]);
        assertEquals("SIGNATURE ETIOLOGY ALLOCATION PERCENT", signAllocTable[1]);
        assertEquals("Sig7 Ultraviolet light exposure 23622 60%", signAllocTable[2]);
        assertEquals("MISALLOC - 4983 13%", signAllocTable[11]);

        // Check the image
        assertEquals(1, page3Ripper.numberOfImages());
        final RectangleInPage topLeftCorner = new RectangleInPage(0.0, 0.0, 0.2, 0.09);
        page3Ripper.checkHasImageWithinBoundsOfGivenSize(topLeftCorner, 1234, 1200);
    }

    private void checkPage4(final PDDocument document) throws IOException
    {
        PageRipper page4Ripper = new PageRipper(document.getPage(3));
        checkHeaderAndFooter(page4Ripper, 4, true);

        // Check the images
        assertEquals(4, page4Ripper.numberOfImages());
        final RectangleInPage topLeftCorner = new RectangleInPage(0.0, 0.0, 0.2, 0.09);
        page4Ripper.checkHasImageWithinBoundsOfGivenSize(topLeftCorner, 1234, 1200);

        final RectangleInPage driver1PlotLocation = new RectangleInPage(0.0, 0.1, 0.5, 0.45);
        page4Ripper.checkHasImageWithinBoundsOfGivenSizeAndColor(driver1PlotLocation, 3000, 3150, new Color(50, 205, 50));

        final RectangleInPage driver2PlotLocation = new RectangleInPage(0.5, 0.1, 1.0, 0.45);
        page4Ripper.checkHasImageWithinBoundsOfGivenSizeAndColor(driver2PlotLocation, 3000, 3150, new Color(0, 191, 255));

        final RectangleInPage driver3PlotLocation = new RectangleInPage(0.0, 0.4, 0.5, 0.8);
        page4Ripper.checkHasImageWithinBoundsOfGivenSizeAndColor(driver3PlotLocation, 3000, 3150, new Color(127, 255, 0));
    }

    private void checkPage5(final PDDocument document) throws IOException
    {
        PageRipper page5Ripper = new PageRipper(document.getPage(4));
        checkHeaderAndFooter(page5Ripper, 5, true);

        String[] wholePage = page5Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.1), new PositionInPage(1.0, 0.6));
        assertEquals(12, wholePage.length);
        assertEquals("Germline Findings", wholePage[0]);
        assertEquals("Small Variants (0)", wholePage[1]);
        assertEquals("NONE", wholePage[2]);
        assertEquals("Amplifications and Deletions (0)", wholePage[3]);
        assertEquals("NONE", wholePage[4]);
        assertEquals("Disruptions (0)", wholePage[5]);
        assertEquals("NONE", wholePage[6]);
        //        assertEquals("Chromosomal Aberrations (0)", wholePage[7]);
        assertEquals("NONE", wholePage[8]);
        assertEquals("Pharmacogenetics (1)", wholePage[9]);
        assertEquals("GENE HAPLOTYPE GENOTYPE FUNCTION LINKED DRUGS SOURCE", wholePage[10]);
        assertEquals("DPYD *1 HOM Normal Function 5-Fluorouracil;Capecitabine;Tegafur PHARMGKB", wholePage[11]);

        // Check the image
        assertEquals(1, page5Ripper.numberOfImages());
        final RectangleInPage topLeftCorner = new RectangleInPage(0.0, 0.0, 0.2, 0.09);
        page5Ripper.checkHasImageWithinBoundsOfGivenSize(topLeftCorner, 1234, 1200);
    }

    private void checkPage6(final PDDocument document) throws IOException
    {
        PageRipper page6Ripper = new PageRipper(document.getPage(5));
        checkHeaderAndFooter(page6Ripper, 6, true);

        String[] immunologyHeading = page6Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.1), new PositionInPage(1.0, 0.145));
        assertEquals(1, immunologyHeading.length);
        assertEquals("Immunology", immunologyHeading[0]);

        String[] classITable = page6Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.145), new PositionInPage(1.0, 0.31));
        assertEquals(8, classITable.length);
        assertEquals("HLA Class I Alleles", classITable[0]);
        assertEquals("ALLELE QC STATUS REF FRAGS TUMOR FRAGS TUMOR CN SOMATIC MUTATIONS", classITable[1]);
        assertEquals("A*01:01 PASS 258 706 2.0 NONE", classITable[2]);
        assertEquals("A*01:01 PASS 258 707 1.4 NONE", classITable[3]);
        assertEquals("B*08:01 PASS 216 640 1.4 NONE", classITable[4]);
        assertEquals("B*40:02 PASS 212 751 2.0 NONE", classITable[5]);
        assertEquals("C*03:04 PASS 249 803 2.0 NONE", classITable[6]);
        assertEquals("C*07:01 PASS 225 621 1.4 NONE", classITable[7]);

        String[] classIITable = page6Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.325), new PositionInPage(1.0, 0.52));
        assertEquals(10, classIITable.length);
        assertEquals("HLA Class II Alleles", classIITable[0]);
        assertEquals("ALLELE QC STATUS REF FRAGS TUMOR FRAGS TUMOR CN SOMATIC MUTATIONS", classIITable[1]);
        assertEquals("DPA1*01:03 PASS 127 440 2.0 NONE", classIITable[2]);
        assertEquals("DPA1*02:01 PASS 103 296 1.4 NONE", classIITable[3]);
        assertEquals("DPB1*01:01 PASS 150 374 1.4 NONE", classIITable[4]);
        assertEquals("DPB1*02:01 PASS 129 491 2.0 NONE", classIITable[5]);
        assertEquals("DQA1*05:01 PASS 92 293 2.0 NONE", classIITable[6]);
        assertEquals("DQA1*05:05 PASS 99 298 1.4 NONE", classIITable[7]);
        assertEquals("DQB1*02:01 PASS 106 314 1.4 NONE", classIITable[8]);
        assertEquals("DQB1*03:01 PASS 113 415 2.0 NONE", classIITable[9]);

        assertEquals(1, page6Ripper.numberOfImages());
        final RectangleInPage topLeftCorner = new RectangleInPage(0.0, 0.0, 0.2, 0.09);
        page6Ripper.checkHasImageWithinBoundsOfGivenSize(topLeftCorner, 1234, 1200);
    }

    private void checkPage7(final PDDocument document) throws IOException
    {
        PageRipper pageRipper = new PageRipper(document.getPage(6));
        checkHeaderAndFooter(pageRipper, 7, false);

        String[] wholePage = pageRipper.getLinesInRectangle(new PositionInPage(0.0, 0.12), new PositionInPage(1.0, 0.2));
        assertEquals(1, wholePage.length);
        assertEquals("Tissue of Origin", wholePage[0]);

        // Check the images - note that this page has landscape orientation
        assertEquals(2, pageRipper.numberOfImages());
        final RectangleInPage topLeftCorner = new RectangleInPage(0.0, 0.0, 0.2, 0.15);
        pageRipper.checkHasImageWithinBoundsOfGivenSize(topLeftCorner, 1234, 1200);

        final RectangleInPage driver1PlotLocation = new RectangleInPage(0.05, 0.1, 0.95, 0.95);
        pageRipper.checkHasImageWithinBoundsOfGivenSizeAndColor(driver1PlotLocation, 6000, 3300, new Color(255, 20, 147));
    }

    private void checkPage8(final PDDocument document) throws IOException
    {
        PageRipper pageRipper = new PageRipper(document.getPage(7));
        checkHeaderAndFooter(pageRipper, 8, false);

        String[] wholePage = pageRipper.getLinesInRectangle(new PositionInPage(0.0, 0.12), new PositionInPage(1.0, 0.2));
        assertEquals(1, wholePage.length);
        assertEquals("Quality Control", wholePage[0]);

        // Check the images - note that this page has landscape orientation
        assertEquals(2, pageRipper.numberOfImages());
        final RectangleInPage topLeftCorner = new RectangleInPage(0.0, 0.0, 0.2, 0.15);
        pageRipper.checkHasImageWithinBoundsOfGivenSize(topLeftCorner, 1234, 1200);

        final RectangleInPage qcImage = new RectangleInPage(0.05, 0.1, 0.95, 0.95);
        pageRipper.checkHasImageWithinBoundsOfGivenSizeAndColor(qcImage, 5100, 3000, new Color(139, 69, 19));
    }

    private void checkPage9(final PDDocument document) throws IOException
    {
        PageRipper page9Ripper = new PageRipper(document.getPage(8));
        checkHeaderAndFooter(page9Ripper, 9, false);

        String[] wholePage = page9Ripper.getLinesInRectangle(new PositionInPage(0.0, 0.12), new PositionInPage(1.0, 0.2));
        assertEquals(1, wholePage.length);
        assertEquals("Purity and Ploidy", wholePage[0]);

        // Check the images - note that this page has landscape orientation
        assertEquals(7, page9Ripper.numberOfImages());
        final RectangleInPage topLeftCorner = new RectangleInPage(0.0, 0.0, 0.2, 0.15);
        page9Ripper.checkHasImageWithinBoundsOfGivenSize(topLeftCorner, 1234, 1200);

        final RectangleInPage pinkPic = new RectangleInPage(0.0, 0.05, 0.4, 0.7);
        page9Ripper.checkHasImageWithinBoundsOfGivenSizeAndColor(pinkPic, 3000, 3000, new Color(255, 105, 180));

        final RectangleInPage redPic = new RectangleInPage(0.2, 0.05, 0.6, 0.7);
        page9Ripper.checkHasImageWithinBoundsOfGivenSizeAndColor(redPic, 1440, 1200, new Color(255, 69, 0));

        final RectangleInPage bluePic = new RectangleInPage(0.5, 0.1, 0.95, 0.7);
        page9Ripper.checkHasImageWithinBoundsOfGivenSizeAndColor(bluePic, 2400, 1800, new Color(30, 144, 255));

        final RectangleInPage tealPic = new RectangleInPage(0.0, 0.5, 0.4, 0.9);
        page9Ripper.checkHasImageWithinBoundsOfGivenSizeAndColor(tealPic, 1440, 1200, new Color(0, 128, 128));

        final RectangleInPage orangePic = new RectangleInPage(0.2, 0.5, 0.6, 0.9);
        page9Ripper.checkHasImageWithinBoundsOfGivenSizeAndColor(orangePic, 1440, 1200, new Color(255, 140, 0));

        final RectangleInPage yellowPic = new RectangleInPage(0.5, 0.5, 0.95, 0.9);
        page9Ripper.checkHasImageWithinBoundsOfGivenSizeAndColor(yellowPic, 2400, 1200, new Color(255, 215, 0));
    }

    private void checkHeaderAndFooter(PageRipper ripper, int page, boolean isPortrait) throws IOException
    {
        double topLeftXMax = isPortrait ? 0.48 : 0.45;
        double topLeftYMax = isPortrait ? 0.05 : 0.12;
        String[] topLeft = ripper.getLinesInRectangle(new PositionInPage(0.0, 0.0), new PositionInPage(topLeftXMax, topLeftYMax));
        assertEquals(1, topLeft.length);
        assertEquals("ORANGE Report", topLeft[0]);

        double topRightXMin = isPortrait ? 0.52 : 0.7;
        double topRightYMax = isPortrait ? 0.05 : 0.12;
        String[] topRight = ripper.getLinesInRectangle(new PositionInPage(topRightXMin, 0.0), new PositionInPage(1.0, topRightYMax));
        assertEquals(2, topRight.length);
        assertEquals("SAMPLE", topRight[0]);
        assertEquals(tumor, topRight[1]);

        String[] bottomLeft = ripper.getLinesInRectangle(new PositionInPage(0.0, 0.95), new PositionInPage(0.2, 1.0));
        assertEquals(1, bottomLeft.length);
        assertEquals(page + "/9", bottomLeft[0]);
    }

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
