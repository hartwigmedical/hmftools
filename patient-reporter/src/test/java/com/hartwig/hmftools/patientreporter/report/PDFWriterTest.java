package com.hartwig.hmftools.patientreporter.report;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.HmfReporterDataLoader;
import com.hartwig.hmftools.patientreporter.ImmutablePatientReport;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableReason;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableStudy;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReportType;
import com.hartwig.hmftools.patientreporter.copynumber.ImmutableCopyNumberReport;
import com.hartwig.hmftools.patientreporter.variants.ImmutableVariantReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.junit.Test;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriterTest {

    private static final boolean SHOW_AND_PRINT = false;
    private static final boolean WRITE_TO_PDF = false;

    private static final String REPORT_BASE_DIR = System.getProperty("user.home");

    @Test
    public void canGeneratePatientReport() throws DRException, IOException, HartwigException {
        final String sample = "CPCT11111111T";
        final FittedPurity fittedPurity =
                ImmutableFittedPurity.builder().purity(0.58).diploidProportion(0).normFactor(0).score(0).ploidy(2).build();
        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, fittedPurity);

        final VariantReport variant1 = ImmutableVariantReport.builder()
                .gene("BRAF")
                .chromosome("7")
                .position(140453136)
                .ref("A")
                .alt("T")
                .transcript("ENST00000377970.6")
                .hgvsCoding("c.1799T>A")
                .hgvsProtein("p.Val600Glu")
                .consequence("missense variant")
                .cosmicID("COSM476")
                .alleleReadCount(18)
                .totalReadCount(99)
                .baf("AAAB")
                .impliedVAF(purityAdjuster.purityAdjustedVAF(4, 0.18 / 0.99))
                .build();
        final VariantReport variant2 = ImmutableVariantReport.builder()
                .gene("MYC")
                .chromosome("8")
                .position(128748854)
                .ref("GG")
                .alt("CA")
                .transcript("ENST00000377970.2")
                .hgvsCoding("c.15_16delinsCA")
                .hgvsProtein("p.Val6Ile")
                .consequence("missense variant")
                .cosmicID("")
                .alleleReadCount(20)
                .totalReadCount(88)
                .impliedVAF(purityAdjuster.purityAdjustedVAF(2, 0.2 / 0.88))
                .baf("AB")
                .build();
        final VariantReport variant3 = ImmutableVariantReport.builder()
                .gene("TP53")
                .chromosome("17")
                .position(7577111)
                .ref("GCACAAA")
                .alt("G")
                .transcript("ENST00000269305.4")
                .hgvsCoding("c.821_826delTTTGTG")
                .hgvsProtein("p.Val274_Cys275del")
                .consequence("inframe deletion")
                .alleleReadCount(20)
                .totalReadCount(87)
                .impliedVAF(purityAdjuster.purityAdjustedVAF(3, 0.20 / 0.87))
                .baf("AAA")
                .build();
        final List<VariantReport> variants = Lists.newArrayList(variant1, variant2, variant3);

        final CopyNumberReport copyNumber1 = ImmutableCopyNumberReport.builder()
                .chromosome("2")
                .chromosomeBand("p23.1-p23.2")
                .gene("ALK")
                .copyNumber(0)
                .type(CopyNumberReportType.LOSS)
                .build();
        final CopyNumberReport copyNumber2 = ImmutableCopyNumberReport.builder()
                .chromosome("3")
                .chromosomeBand("q26.32")
                .gene("PIK3CA")
                .copyNumber(9)
                .type(CopyNumberReportType.GAIN)
                .build();
        final List<CopyNumberReport> copyNumbers = Lists.newArrayList(copyNumber1, copyNumber2);

        final int mutationalLoad = 361;
        final String tumorType = "Melanoma";
        final Double pathologyTumorPercentage = 0.6;

        final PatientReport patientReport =
                ImmutablePatientReport.of(sample, variants, copyNumbers, mutationalLoad, tumorType, pathologyTumorPercentage, "58%",
                        "FC000001", "CSB000000", LocalDate.parse("2016-01-05"), LocalDate.parse("2016-01-01"));

        final String genePanelPath = Resources.getResource("bed").getPath() + File.separator + "hmf_gene_panel.tsv";
        final String drupFilterPath = Resources.getResource("csv").getPath() + File.separator + "drup_genes.csv";
        final String cosmicPath = Resources.getResource("csv").getPath() + File.separator + "cosmic_slice.csv";
        final String centerPath = Resources.getResource("center").getPath() + File.separator + "centers.csv";
        final String signaturePath = Resources.getResource("signature").getPath() + File.separator + "signature.png";

        // KODU: Refers to the actual cosmic path on datastore:
        // final String cosmicPath = REPORT_BASE_DIR + "/hmf/tmp/170529_grch37_cosmic_census.csv";
        final HmfReporterData reporterData =
                HmfReporterDataLoader.buildFromFiles(genePanelPath, drupFilterPath, cosmicPath, centerPath, signaturePath);

        final InputStream logoStream = Resources.asByteSource(Resources.getResource(PDFWriter.REPORT_LOGO_PATH)).openStream();
        final JasperReportBuilder report = PDFWriter.generatePatientReport(patientReport, logoStream, reporterData);
        assertNotNull(report);

        if (SHOW_AND_PRINT) {
            report.show().print();
        }

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + "/hmf/tmp/test_report.pdf"));
        }
        logoStream.close();
    }

    @Test
    public void canGenerateNotSequenceableReport() throws DRException, IOException {
        final String sample = "CPCT11111111T";
        final String tumorType = "Melanoma";
        final NotSequenceableReason reason = NotSequenceableReason.LOW_TUMOR_PERCENTAGE;
        final String tumorPercentageString = "10%";
        final InputStream logoStream = Resources.asByteSource(Resources.getResource(PDFWriter.REPORT_LOGO_PATH)).openStream();
        final NotSequenceableStudy study = NotSequenceableStudy.CPCT;

        final JasperReportBuilder report =
                PDFWriter.generateNotSequenceableReport(sample, tumorType, tumorPercentageString, reason, study, logoStream);
        assertNotNull(report);

        if (SHOW_AND_PRINT) {
            report.show().print();
        }

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + "/hmf/tmp/low_tumor_percentage_report.pdf"));
        }
        logoStream.close();
    }
}