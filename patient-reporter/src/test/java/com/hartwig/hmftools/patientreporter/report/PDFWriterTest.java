package com.hartwig.hmftools.patientreporter.report;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableReason;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriterTest {

    private static final boolean SHOW_AND_PRINT = false;
    private static final boolean WRITE_TO_PDF = false;

    private static final String RESOURCE_PATH = Resources.getResource("pdf").getPath();
    private static final String REPORT_LOGO = RESOURCE_PATH + File.separator + "hartwig_logo.jpg";

    @Test
    public void canGeneratePatientReport() throws DRException, IOException, EmptyFileException {
        final String sample = "CPCT11111111T";
        final VariantReport variant1 = new VariantReport.Builder().gene("BRAF").position("7:140453136").ref("A").alt(
                "T").transcript("ENST00000288602.6").hgvsCoding("c.1799T>A").hgvsProtein("p.Val600Glu").consequence(
                "missense variant").cosmicID("COSM476").alleleReadCount(34).totalReadCount(99).build();
        final VariantReport variant2 = new VariantReport.Builder().gene("FGFR3").position("4:1795557").ref("T").alt(
                "C").transcript("ENST00000340107.4").hgvsCoding("c.-102-3T>C").hgvsProtein("").consequence(
                "splice region variant").cosmicID("").alleleReadCount(12).totalReadCount(88).build();
        final VariantReport variant3 = new VariantReport.Builder().gene("GNAQ").position("9:80409488").ref("T").alt(
                "G").transcript("ENST00000286548.4").hgvsCoding("c.626A>C").hgvsProtein("p.Gln209Pro").consequence(
                "missense variant").cosmicID("COSM28758").alleleReadCount(57).totalReadCount(121).build();
        final List<VariantReport> variants = Lists.newArrayList(variant1, variant2, variant3);

        final CopyNumberReport copyNumber1 = new CopyNumberReport.Builder().chromosome("3").gene("PIK3CA").transcript(
                "ENST00000263967.3").copyNumber(6).build();
        final CopyNumberReport copyNumber2 = new CopyNumberReport.Builder().chromosome("2").gene("ALK").transcript(
                "ENST00000389048.3").copyNumber(0).build();
        final List<CopyNumberReport> copyNumbers = Lists.newArrayList(copyNumber1, copyNumber2);

        final int mutationalLoad = 361;
        final String tumorType = "Melanoma";
        final double tumorPercentage = 0.6;

        final PatientReport patientReport = new PatientReport(sample, variants, copyNumbers, mutationalLoad, tumorType,
                tumorPercentage);

        final JasperReportBuilder report = PDFWriter.generatePatientReport(patientReport, REPORT_LOGO,
                createHMFSlicingRegion());
        assertNotNull(report);

        if (SHOW_AND_PRINT) {
            report.show().print();
        }

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream("/Users/kduyvesteyn/hmf/tmp/test_report.pdf"));
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
            report.toPdf(new FileOutputStream("/Users/kduyvesteyn/hmf/tmp/low_tumor_percentage_report.pdf"));
        }
    }

    @NotNull
    private static Slicer createHMFSlicingRegion() throws IOException, EmptyFileException {
        final String resourcePath = Resources.getResource("bed").getPath();
        return SlicerFactory.fromBedFile(resourcePath + File.separator + "HMF_Slicing.bed");
    }
}