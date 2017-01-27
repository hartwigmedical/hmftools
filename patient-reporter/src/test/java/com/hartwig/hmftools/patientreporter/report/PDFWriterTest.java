package com.hartwig.hmftools.patientreporter.report;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.junit.Test;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriterTest {

    private static final String RESOURCE_PATH = Resources.getResource("pdf").getPath();
    private static final String HMF_LOGO = RESOURCE_PATH + File.separator + "hartwig_logo.jpg";

    @Test
    public void canGenerateReport() throws DRException, FileNotFoundException {
        final String sample = "CPCT11111111T";
        final VariantReport variant = new VariantReport.Builder().gene("A").position("3:178936082").ref("G").alt(
                "A").transcript("ENST00000263967.3").hgvsCoding("c.1624G>A").hgvsProtein("p.Glu542Lys").consequence(
                "missense variant").cosmicID("COSM125369").alleleReadCount(75).totalReadCount(151).build();
        final List<VariantReport> variants = Lists.newArrayList(variant, variant, variant);

        final CopyNumberReport copyNumber = new CopyNumberReport.Builder().gene("PIK3CA").transcript(
                "ENST00000263967.3").copyNumber(4).build();
        final List<CopyNumberReport> copyNumbers = Lists.newArrayList(copyNumber, copyNumber);

        final int mutationalLoad = 361;
        final String tumorType = "Melanoma";

        final PatientReport report = new PatientReport(sample, variants, copyNumbers, mutationalLoad, tumorType);

        final JasperReportBuilder pdf = PDFWriter.generatePatientReport(report, HMF_LOGO);
        assertNotNull(pdf);

        // KODU: If you want to visually inspect the report, uncomment the below line!
//        pdf.show().print();
    }
}