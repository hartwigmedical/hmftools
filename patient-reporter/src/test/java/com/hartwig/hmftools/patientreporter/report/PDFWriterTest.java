package com.hartwig.hmftools.patientreporter.report;

import static com.hartwig.hmftools.common.purity.PurityAdjustment.purityAdjustedVAF;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableReason;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;
import com.hartwig.hmftools.patientreporter.genePanel.GenePanelModel;
import com.hartwig.hmftools.patientreporter.variants.ImmutableVariantReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriterTest {

    private static final boolean SHOW_AND_PRINT = false;
    private static final boolean WRITE_TO_PDF = false;

    private static final String RESOURCE_PATH = Resources.getResource("pdf").getPath();
    private static final String CSV_PATH = Resources.getResource("csv").getPath();
    private static final String REPORT_LOGO = RESOURCE_PATH + File.separator + "hartwig_logo.jpg";
    private static final String DRUP_GENES_CSV = RESOURCE_PATH + File.separator + "drup_genes.csv";
    private static final String GENE_PANEL_CSV = CSV_PATH + File.separator + "gene_panel_data.csv";

    @Test
    public void canGeneratePatientReport() throws DRException, IOException, HartwigException {
        final String sample = "CPCT11111111T";
        final FittedPurity fittedPurity = ImmutableFittedPurity.builder().purity(0.58).diploidProportion(
                0).modelBAFDeviation(0).normFactor(0).score(0).build();

        final VariantReport variant1 = ImmutableVariantReport.builder().gene("BRAF").chromosome("7").position(
                140453136).ref("A").alt("T").transcript("ENST00000377970.6").hgvsCoding("c.1799T>A").hgvsProtein(
                "p.Val600Glu").consequence("missense variant").cosmicID("COSM476").alleleReadCount(34).totalReadCount(
                99).impliedVAF(purityAdjustedVAF(fittedPurity.purity(), 4, 0.34 / 0.99)).baf("AAAB").build();
        final VariantReport variant2 = ImmutableVariantReport.builder().gene("MYC").chromosome("8").position(
                128748854).ref("GG").alt("CA").transcript("ENST00000377970.2").hgvsCoding(
                "c.15_16delGGinsCA").hgvsProtein("p.ArgVal5ArgIle").consequence("missense variant").cosmicID(
                "").alleleReadCount(12).totalReadCount(88).impliedVAF(
                purityAdjustedVAF(fittedPurity.purity(), 2, 0.12 / 0.88)).baf("AB").build();
        final VariantReport variant3 = ImmutableVariantReport.builder().gene("TP53").chromosome("17").position(
                7577111).ref("GCACAAA").alt("G").transcript("ENST00000269305.4").hgvsCoding(
                "c.821_826delTTTGTG").hgvsProtein("p.Val274_Cys275del").consequence(
                "inframe deletion").alleleReadCount(21).totalReadCount(87).impliedVAF(
                purityAdjustedVAF(fittedPurity.purity(), 3, 0.21 / 0.87)).baf("AAA").build();
        final List<VariantReport> variants = Lists.newArrayList(variant1, variant2, variant3);

        final CopyNumberReport copyNumber1 = new CopyNumberReport.Builder().chromosome("2").gene("ALK").transcript(
                "ENST00000389048.3").copyNumber(0).build();
        final CopyNumberReport copyNumber2 = new CopyNumberReport.Builder().chromosome("3").gene("PIK3CA").transcript(
                "ENST00000263967.3").copyNumber(6).build();
        final List<CopyNumberReport> copyNumbers = Lists.newArrayList(copyNumber1, copyNumber2);

        final int mutationalLoad = 361;
        final String tumorType = "Melanoma";
        final Double pathologyTumorPercentage = 0.6;

        final PatientReport patientReport = new PatientReport(sample, variants, copyNumbers, mutationalLoad, tumorType,
                pathologyTumorPercentage, fittedPurity);
        final DrupFilter drupFilter = new DrupFilter(DRUP_GENES_CSV);
        final GenePanelModel genePanelModel = new GenePanelModel(GENE_PANEL_CSV);

        final JasperReportBuilder report = PDFWriter.generatePatientReport(patientReport, REPORT_LOGO,
                createHMFSlicingRegion(), drupFilter, genePanelModel);
        assertNotNull(report);

        if (SHOW_AND_PRINT) {
            report.show().print();
        }

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream("/Users/korneel/hmf/tmp/test_report.pdf"));
        }
    }

    @Test
    public void canGenerateNotSequenceableReport() throws DRException, FileNotFoundException {
        final String sample = "CPCT11111111T";
        final String tumorType = "Melanoma";
        final NotSequenceableReason reason = NotSequenceableReason.LOW_TUMOR_PERCENTAGE;
        final String tumorPercentageString = "10%";

        final JasperReportBuilder report = PDFWriter.generateNotSequenceableReport(sample, tumorType,
                tumorPercentageString, reason, REPORT_LOGO);
        assertNotNull(report);

        if (SHOW_AND_PRINT) {
            report.show().print();
        }

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream("/Users/korneel/hmf/tmp/low_tumor_percentage_report.pdf"));
        }
    }

    @NotNull
    private static Slicer createHMFSlicingRegion() throws IOException, EmptyFileException {
        final String resourcePath = Resources.getResource("bed").getPath();
        return SlicerFactory.fromBedFile(resourcePath + File.separator + "HMF_Slicing_v2.bed");
    }
}