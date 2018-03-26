package com.hartwig.hmftools.patientreporter.report;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.mockedAlterations;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testBaseReporterData;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testHmfReporterData;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Locale;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.patientreporter.BaseReporterData;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.ImmutableNotSequencedPatientReport;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.ImmutableSequencedPatientReport;
import com.hartwig.hmftools.patientreporter.NotSequencedPatientReport;
import com.hartwig.hmftools.patientreporter.PatientReporterTestUtil;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.SequencedPatientReport;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableReason;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableStudy;
import com.hartwig.hmftools.patientreporter.report.data.Alteration;
import com.hartwig.hmftools.patientreporter.report.data.GeneDisruptionData;
import com.hartwig.hmftools.patientreporter.report.data.GeneFusionData;
import com.hartwig.hmftools.patientreporter.report.data.ImmutableGeneDisruptionData;
import com.hartwig.hmftools.patientreporter.report.data.ImmutableGeneFusionData;
import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;
import com.hartwig.hmftools.patientreporter.variants.ImmutableVariantReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriterTest {

    private static final boolean RUN_CIVIC_ANALYSIS = false;
    private static final boolean SHOW_AND_PRINT = false;
    private static final boolean WRITE_TO_PDF = false;

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";

    private static final DateTimeFormatter FORMATTER = DateTimeFormatter.ofPattern("dd-MMM-yyyy", Locale.ENGLISH);

    @Test
    public void canGenerateSequenceReport() throws DRException, IOException {
        final double pathologyTumorPercentage = 0.6;
        final double impliedTumorPurity = 0.58;
        final int mutationalLoad = 361;
        final double microsatelliteIndicator = 2.1;

        final HmfReporterData reporterData = testHmfReporterData();
        final BaseReporterData baseReporterData = testBaseReporterData();
        final FittedPurity fittedPurity = createFittedPurity(impliedTumorPurity);

        final List<VariantReport> variants = createTestVariants(new PurityAdjuster(Gender.MALE, fittedPurity));
        final List<GeneCopyNumber> copyNumbers = createTestCopyNumbers();
        final List<GeneDisruptionData> disruptions = createTestDisruptions();
        final List<GeneFusionData> fusions = createTestFusions();

        final SampleReport sampleReport = testSampleReport(pathologyTumorPercentage);
        final List<Alteration> alterations = RUN_CIVIC_ANALYSIS ? PatientReporterTestUtil.runCivicAnalysis(variants,
                copyNumbers,
                disruptions,
                fusions, reporterData.panelGeneModel(), sampleReport.cancerType()) : mockedAlterations();

        final SequencedPatientReport patientReport = ImmutableSequencedPatientReport.of(sampleReport,
                variants,
                mutationalLoad,
                microsatelliteIndicator,
                copyNumbers,
                disruptions,
                fusions,
                PatientReportFormat.formatPercent(impliedTumorPurity),
                alterations,
                Resources.getResource("circos" + File.separator + "circos_example.png").getPath(),
                Optional.of("this is a test report and does not relate to any real CPCT patient"),
                baseReporterData.signaturePath());

        final JasperReportBuilder mainReport = PDFWriter.generatePatientReport(patientReport, reporterData);
        assertNotNull(mainReport);

        final JasperReportBuilder evidenceReport = EvidenceReport.generate(patientReport);
        assertNotNull(evidenceReport);

        if (SHOW_AND_PRINT) {
            mainReport.show().print();
        }

        if (WRITE_TO_PDF) {
            mainReport.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "test_report.pdf"));
            evidenceReport.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "test_evidence_report.pdf"));
        }
    }

    @NotNull
    private static FittedPurity createFittedPurity(double impliedPurity) {
        return ImmutableFittedPurity.builder().purity(impliedPurity).diploidProportion(0).normFactor(0).score(0).ploidy(2).build();
    }

    @NotNull
    private static List<VariantReport> createTestVariants(@NotNull final PurityAdjuster purityAdjuster) {
        final VariantReport variant1 = ImmutableVariantReport.builder()
                .gene("BRAF")
                .variant(createTestVariant("7", 140453136, "A", "T"))
                .transcript("ENST00000377970.6")
                .hgvsCoding("c.1799T>A")
                .hgvsProtein("p.Val600Glu")
                .consequence("missense variant")
                .cosmicID("COSM476")
                .alleleReadCount(18)
                .totalReadCount(99)
                .baf("AAAB")
                .impliedVAF(purityAdjuster.purityAdjustedVAF("7", 4, 0.18 / 0.99))
                .build();

        final VariantReport variant2 = ImmutableVariantReport.builder()
                .gene("MYC")
                .variant(createTestVariant("8", 128748854, "GG", "CA"))
                .transcript("ENST00000377970.2")
                .hgvsCoding("c.15_16delinsCA")
                .hgvsProtein("p.Val6Ile")
                .consequence("missense variant")
                .cosmicID("")
                .alleleReadCount(20)
                .totalReadCount(88)
                .impliedVAF(purityAdjuster.purityAdjustedVAF("8", 2, 0.2 / 0.88))
                .baf("AB")
                .build();

        final VariantReport variant3 = ImmutableVariantReport.builder()
                .gene("TP53")
                .variant(createTestVariant("17", 7577111, "GCACAAA", "G"))
                .transcript("ENST00000269305.4")
                .hgvsCoding("c.821_826delTTTGTG")
                .hgvsProtein("p.Val274_Cys275del")
                .consequence("inframe deletion")
                .alleleReadCount(20)
                .totalReadCount(87)
                .impliedVAF(purityAdjuster.purityAdjustedVAF("17", 3, 0.20 / 0.87))
                .baf("AAA")
                .build();
        return Lists.newArrayList(variant1, variant2, variant3);
    }

    @NotNull
    private static SomaticVariant createTestVariant(@NotNull final String chromosome, final long position, @NotNull final String ref,
            @NotNull final String alt) {
        return SomaticVariantTestBuilderFactory.create().chromosome(chromosome).position(position).ref(ref).alt(alt).build();
    }

    @NotNull
    private static List<GeneCopyNumber> createTestCopyNumbers() {
        final GeneCopyNumber copyNumber1 =
                createCopyNumberBuilder().chromosome("9").chromosomeBand("p21.3").gene("CDKN2A").minCopyNumber(0).build();
        final GeneCopyNumber copyNumber2 =
                createCopyNumberBuilder().chromosome("17").chromosomeBand("q12").gene("ERBB2").minCopyNumber(9).build();
        return Lists.newArrayList(copyNumber1, copyNumber2);
    }

    @NotNull
    private static ImmutableGeneCopyNumber.Builder createCopyNumberBuilder() {
        return ImmutableGeneCopyNumber.builder()
                .start(1)
                .end(2)
                .minRegionStart(0)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEnd(0)
                .minRegionEndSupport(SegmentSupport.NONE)
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegions(1)
                .germlineHet2HomRegions(0)
                .germlineHomRegions(0)
                .somaticRegions(1)
                .maxCopyNumber(0)
                .transcriptID("trans")
                .transcriptVersion(0)
                .nonsenseBiallelicCount(0)
                .nonsenseNonBiallelicCount(0)
                .nonsenseNonBiallelicPloidy(0)
                .spliceBiallelicCount(0)
                .spliceNonBiallelicCount(0)
                .spliceNonBiallelicPloidy(0)
                .missenseBiallelicCount(0)
                .missenseNonBiallelicCount(0)
                .missenseNonBiallelicPloidy(0)
                .minMinorAllelePloidy(0);
    }

    @NotNull
    private static List<GeneFusionData> createTestFusions() {
        return Lists.newArrayList(ImmutableGeneFusionData.builder()
                        .geneStart("TMPRSS2")
                        .geneStartTranscript("ENST00000398585.7")
                        .geneStartEntrezIds(Lists.newArrayList(7113))
                        .geneContextStart("Exon 1")
                        .geneEnd("PNPLA7")
                        .geneEndTranscript("ENST00000406427.5")
                        .geneEndEntrezIds(Lists.newArrayList(375775))
                        .geneContextEnd("Exon 13")
                        .copies("1.0")
                        .cosmicURL("")
                        .build(),
                ImmutableGeneFusionData.builder()
                        .geneStart("CLCN6")
                        .geneStartTranscript("ENST00000346436.10")
                        .geneStartEntrezIds(Lists.newArrayList(1185))
                        .geneContextStart("Exon 1")
                        .geneEnd("BRAF")
                        .geneEndTranscript("ENST00000288602.10")
                        .geneEndEntrezIds(Lists.newArrayList(673))
                        .geneContextEnd("Exon 13")
                        .copies("1.0")
                        .cosmicURL("http://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=2&gid=54500")
                        .build());
    }

    @NotNull
    private static List<GeneDisruptionData> createTestDisruptions() {
        final GeneDisruptionData disruption1 = ImmutableGeneDisruptionData.builder()
                .chromosome("2")
                .gene("ERBB4")
                .geneContext("Intron 4")
                .type("INV")
                .copies("1.0").chromosomeBand("q34")
                .build();

        final GeneDisruptionData disruption2 = ImmutableGeneDisruptionData.builder()
                .chromosome("2")
                .gene("ERBB4")
                .geneContext("Intron 20")
                .type("INV")
                .copies("1.0").chromosomeBand("q34")
                .build();

        final GeneDisruptionData disruption3 = ImmutableGeneDisruptionData.builder()
                .chromosome("3")
                .gene("PIK3CB")
                .geneContext("Intron 1")
                .type("INS")
                .copies("3.0").chromosomeBand("q22.3")
                .build();

        final GeneDisruptionData disruption4 = ImmutableGeneDisruptionData.builder()
                .chromosome("8")
                .gene("NRG1")
                .geneContext("Intron 1")
                .type("DUP")
                .copies("0.3").chromosomeBand("p12")
                .build();

        final GeneDisruptionData disruption5 = ImmutableGeneDisruptionData.builder()
                .chromosome("8")
                .gene("NRG1")
                .geneContext("Intron 1")
                .type("DEL")
                .copies("0.2").chromosomeBand("p12")
                .build();

        final GeneDisruptionData disruption6 = ImmutableGeneDisruptionData.builder()
                .chromosome("17")
                .gene("CDK12")
                .geneContext("Intron 12")
                .type("BND")
                .copies("1.0").chromosomeBand("q12")
                .build();

        return Lists.newArrayList(disruption1, disruption2, disruption3, disruption4, disruption5, disruption6);
    }

    @Test
    public void canGenerateLowTumorPercentageReport() throws DRException, IOException {
        final JasperReportBuilder report = generateNotSequenceableCPCTReport(0.1, NotSequenceableReason.LOW_TUMOR_PERCENTAGE);
        assertNotNull(report);

        if (SHOW_AND_PRINT) {
            report.show().print();
        }

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_low_tumor_percentage_report.pdf"));
        }
    }

    @Test
    public void canGenerateLowDNAYieldReport() throws DRException, IOException {
        final JasperReportBuilder report = generateNotSequenceableCPCTReport(0.6, NotSequenceableReason.LOW_DNA_YIELD);
        assertNotNull(report);

        if (SHOW_AND_PRINT) {
            report.show().print();
        }

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_low_dna_yield_report.pdf"));
        }
    }

    @Test
    public void canGeneratePostDNAIsolationFailReport() throws DRException, IOException {
        final JasperReportBuilder report = generateNotSequenceableCPCTReport(0.6, NotSequenceableReason.POST_ISOLATION_FAIL);
        assertNotNull(report);

        if (SHOW_AND_PRINT) {
            report.show().print();
        }

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_post_dna_isolation_fail_report.pdf"));
        }
    }

    @NotNull
    private static JasperReportBuilder generateNotSequenceableCPCTReport(final double pathologyTumorEstimate,
            @NotNull final NotSequenceableReason reason) throws IOException {
        final NotSequencedPatientReport patientReport = ImmutableNotSequencedPatientReport.of(testSampleReport(pathologyTumorEstimate),
                reason,
                NotSequenceableStudy.CPCT,
                Optional.empty(),
                PatientReporterTestUtil.SIGNATURE_PATH);

        return PDFWriter.generateNotSequenceableReport(patientReport);
    }

    @NotNull
    private static SampleReport testSampleReport(final double pathologyTumorPercentage) throws IOException {
        final String sample = "CPCT02991111T";
        return ImmutableSampleReport.of(sample,
                "Melanoma",
                pathologyTumorPercentage,
                LocalDate.parse("05-Jan-2016", FORMATTER),
                LocalDate.parse("01-Jan-2016", FORMATTER),
                "PREP013V23-QC037V20-SEQ008V25",
                testBaseReporterData().centerModel().getAddresseeStringForSample(sample));
    }
}