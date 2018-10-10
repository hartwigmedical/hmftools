package com.hartwig.hmftools.patientreporter.report;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testBaseReportData;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testSequencedReportData;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.actionability.somaticvariant.EvidenceItem;
import com.hartwig.hmftools.common.actionability.somaticvariant.VariantEvidenceItems;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.OncoDrivers;
import com.hartwig.hmftools.common.drivercatalog.TsgDrivers;
import com.hartwig.hmftools.common.ecrf.projections.ImmutablePatientTumorLocation;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.Clonality;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.BaseReportData;
import com.hartwig.hmftools.patientreporter.ImmutableAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ImmutableNotAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.NotAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.SequencedReportData;
import com.hartwig.hmftools.patientreporter.algo.NotAnalysableReason;
import com.hartwig.hmftools.patientreporter.algo.NotAnalysableStudy;
import com.hartwig.hmftools.svannotation.annotations.GeneAnnotation;
import com.hartwig.hmftools.svannotation.annotations.GeneDisruption;
import com.hartwig.hmftools.svannotation.annotations.GeneFusion;
import com.hartwig.hmftools.svannotation.annotations.ImmutableGeneDisruption;
import com.hartwig.hmftools.svannotation.annotations.ImmutableGeneFusion;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriterTest {

    private static final boolean WRITE_TO_PDF = false;

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";

    private static final DateTimeFormatter FORMATTER = DateTimeFormatter.ofPattern("dd-MMM-yyyy", Locale.ENGLISH);

    @Test
    public void canGenerateSequenceReport() throws DRException, IOException {
        final double pathologyTumorPercentage = 0.6;
        final double impliedTumorPurity = 0.58;
        final int tumorMutationalLoad = 361;
        final double tumorMutationalBurden = 10.1;
        final double microsatelliteIndelsPerMb = 2.1;

        final SequencedReportData reporterData = testSequencedReportData();
        final BaseReportData baseReportData = testBaseReportData();
        final FittedPurity fittedPurity = createFittedPurity(impliedTumorPurity);

        final List<EnrichedSomaticVariant> variants = createTestVariants(new PurityAdjuster(Gender.MALE, fittedPurity));
        final List<GeneCopyNumber> copyNumbers = createTestCopyNumbers();
        final List<GeneFusion> fusions = createTestFusions();
        final List<GeneDisruption> disruptions = createTestDisruptions();

        final List<DriverCatalog> driverCatalog = Lists.newArrayList();
        driverCatalog.addAll(OncoDrivers.drivers(DndsDriverGeneLikelihoodSupplier.oncoLikelihood(), variants));
        driverCatalog.addAll(TsgDrivers.drivers(DndsDriverGeneLikelihoodSupplier.tsgLikelihood(), variants));

        final List<EvidenceItem> actionVariant = Lists.newArrayList();
        final Map<EnrichedSomaticVariant, VariantEvidenceItems> labelEvidence = Maps.newHashMap();

        final SampleReport sampleReport = testSampleReport(pathologyTumorPercentage);

        final AnalysedPatientReport patientReport = ImmutableAnalysedPatientReport.of(sampleReport,
                FittedPurityStatus.NORMAL,
                impliedTumorPurity,
                variants,
                actionVariant,labelEvidence,
                driverCatalog,
                microsatelliteIndelsPerMb,
                tumorMutationalLoad,
                tumorMutationalBurden,
                Lists.newArrayList(),
                copyNumbers,
                fusions,
                disruptions,
                Resources.getResource("circos/circos_example.png").getPath(),
                Optional.of("this is a test report and does not relate to any real CPCT patient"),
                baseReportData.signaturePath());

        final JasperReportBuilder mainReport = PDFWriter.generatePatientReport(patientReport, reporterData);
        assertNotNull(mainReport);

        if (WRITE_TO_PDF) {
            mainReport.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_test_sequence_report.pdf"));
        }
    }

    @NotNull
    private static FittedPurity createFittedPurity(double impliedPurity) {
        return ImmutableFittedPurity.builder()
                .purity(impliedPurity)
                .diploidProportion(0)
                .normFactor(0)
                .score(0)
                .ploidy(2)
                .somaticDeviation(0)
                .build();
    }

    @NotNull
    private static List<EnrichedSomaticVariant> createTestVariants(@NotNull final PurityAdjuster purityAdjuster) {
        final EnrichedSomaticVariant variant1 = createSomaticVariantBuilder().gene("BRAF")
                .canonicalHgvsCodingImpact("c.1799T>A")
                .canonicalHgvsProteinImpact("p.Val600Glu")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .type(VariantType.SNP)
                .hotspot(Hotspot.HOTSPOT)
                .clonality(Clonality.CLONAL)
                .alleleReadCount(18)
                .totalReadCount(99)
                .adjustedCopyNumber(4)
                .minorAllelePloidy(1)
                .adjustedVAF(purityAdjuster.purityAdjustedVAF("7", 4, 0.18 / 0.99))
                .build();

        final EnrichedSomaticVariant variant2 = createSomaticVariantBuilder().gene("TP53")
                .canonicalHgvsCodingImpact("c.821_826delTTTGTG")
                .canonicalHgvsProteinImpact("p.Val274_Cus275del")
                .canonicalCodingEffect(CodingEffect.NONSENSE_OR_FRAMESHIFT)
                .type(VariantType.INDEL)
                .clonality(Clonality.CLONAL)
                .alleleReadCount(20)
                .totalReadCount(87)
                .adjustedCopyNumber(3)
                .minorAllelePloidy(0)
                .adjustedVAF(purityAdjuster.purityAdjustedVAF("17", 3, 0.20 / 0.87))
                .build();

        final EnrichedSomaticVariant variant3 = createSomaticVariantBuilder().gene("MYC")
                .canonicalHgvsCodingImpact("c.15_16delinsCA")
                .canonicalHgvsProteinImpact("p.Val6Ile")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .type(VariantType.MNP)
                .clonality(Clonality.CLONAL)
                .alleleReadCount(20)
                .totalReadCount(88)
                .adjustedCopyNumber(2)
                .minorAllelePloidy(1)
                .adjustedVAF(purityAdjuster.purityAdjustedVAF("8", 2, 0.2 / 0.88))
                .build();

        return Lists.newArrayList(variant1, variant2, variant3);
    }

    @NotNull
    private static List<GeneCopyNumber> createTestCopyNumbers() {
        final GeneCopyNumber copyNumber1 =
                createCopyNumberBuilder().chromosome("9").chromosomeBand("p21.3").gene("CDKN2A").minCopyNumber(0).maxCopyNumber(0).build();
        final GeneCopyNumber copyNumber2 =
                createCopyNumberBuilder().chromosome("17").chromosomeBand("p13.1").gene("TP53").minCopyNumber(0).maxCopyNumber(2).build();
        final GeneCopyNumber copyNumber3 =
                createCopyNumberBuilder().chromosome("17").chromosomeBand("q12").gene("ERBB2").minCopyNumber(11).maxCopyNumber(11).build();
        return Lists.newArrayList(copyNumber1, copyNumber2, copyNumber3);
    }

    @NotNull
    private static List<GeneFusion> createTestFusions() {
        GeneFusion fusion1 =
                createFusion("TMPRSS2", "ENST00000398585", 4, 5, "PNPLA7", "ENST00000406427", 2, 3, KnownFusionsModel.CIVIC, 0.4);
        GeneFusion fusion2 = createFusion("CLCN6", "ENST00000346436", 1, 2, "BRAF", "ENST00000288602", 8, 9, KnownFusionsModel.ONCOKB, 1D);

        return Lists.newArrayList(fusion1, fusion2);
    }

    @NotNull
    private static GeneFusion createFusion(@NotNull String startGene, @NotNull String startTranscript, int startExonUpstream,
            int startExonDownstream, @NotNull String endGene, @NotNull String endTranscript, int endExonUpstream, int endExonDownstream,
            @NotNull String source, double ploidy) {
        Transcript upstreamTranscript = createFusionLeg(true, startGene, startTranscript, startExonUpstream, startExonDownstream, ploidy);
        Transcript downstreamTranscript = createFusionLeg(false, endGene, endTranscript, endExonUpstream, endExonDownstream, ploidy);

        return ImmutableGeneFusion.builder()
                .reportable(true)
                .primarySource(source)
                .upstreamLinkedAnnotation(upstreamTranscript)
                .downstreamLinkedAnnotation(downstreamTranscript)
                .build();
    }

    @NotNull
    private static Transcript createFusionLeg(boolean isUpstream, @NotNull String gene, @NotNull String transcript, int exonUpstream,
            int exonDownstream, double ploidy) {
        EnrichedStructuralVariant variant = createEnrichedStructuralVariantBuilder().type(StructuralVariantType.BND)
                .start(createEnrichedStructuralVariantLegBuilder().chromosome("any").build())
                .ploidy(ploidy)
                .qualityScore(0)
                .build();
        GeneAnnotation upstreamGene =
                new GeneAnnotation(variant, isUpstream, gene, Strings.EMPTY, 1, Lists.newArrayList(), Lists.newArrayList(), Strings.EMPTY);
        return new Transcript(upstreamGene, transcript, exonUpstream, -1, exonDownstream, -1, 10, true, null, null);
    }

    @NotNull
    private static List<GeneDisruption> createTestDisruptions() {
        GeneDisruption disruption1 = createDisruption(StructuralVariantType.DUP, "8", "p12", "NRG1", 1, 2, 0.3);
        GeneDisruption disruption2 = createDisruption(StructuralVariantType.INV, "2", "q34", "ERBB4", 4, 5, 1D);
        GeneDisruption disruption3 = createDisruption(StructuralVariantType.INS, "3", "q22.3", "PIK3CB", 1, 2, 3D);
        GeneDisruption disruption4 = createDisruption(StructuralVariantType.DEL, "8", "p12", "NRG1", 1, 2, 0.2);
        GeneDisruption disruption5 = createDisruption(StructuralVariantType.BND, "17", "q12", "CDK12", 12, 13, 1D);
        GeneDisruption disruption6 = createDisruption(StructuralVariantType.INV, "2", "q34", "ERBB4", 20, 21, 1D);

        return Lists.newArrayList(disruption1, disruption2, disruption3, disruption4, disruption5, disruption6);
    }

    @NotNull
    private static GeneDisruption createDisruption(@NotNull StructuralVariantType type, @NotNull String chromosome,
            @NotNull String chromosomeBand, @NotNull String gene, int exonUpstream, int exonDownstream, double ploidy) {
        EnrichedStructuralVariantLeg start = createEnrichedStructuralVariantLegBuilder().chromosome(chromosome).build();
        EnrichedStructuralVariant variant = createEnrichedStructuralVariantBuilder().type(type).start(start).ploidy(ploidy).qualityScore(0).build();

        GeneAnnotation geneAnnotation =
                new GeneAnnotation(variant, true, gene, "id", 1, Lists.newArrayList(), Lists.newArrayList(), chromosomeBand);

        Transcript transcript = new Transcript(geneAnnotation, "trans", exonUpstream, -1, exonDownstream, -1, 5, true, null, null);

        return ImmutableGeneDisruption.builder().reportable(true).linkedAnnotation(transcript).build();
    }

    @Test
    public void canGenerateLowTumorPercentageReport() throws DRException, IOException {
        final JasperReportBuilder report = generateNotAnalysableCPCTReport(0.1, NotAnalysableReason.LOW_TUMOR_PERCENTAGE);
        assertNotNull(report);

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_low_tumor_percentage_report.pdf"));
        }
    }

    @Test
    public void canGenerateLowDNAYieldReport() throws DRException, IOException {
        final JasperReportBuilder report = generateNotAnalysableCPCTReport(0.6, NotAnalysableReason.LOW_DNA_YIELD);
        assertNotNull(report);

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_low_dna_yield_report.pdf"));
        }
    }

    @Test
    public void canGeneratePostDNAIsolationFailReport() throws DRException, IOException {
        final JasperReportBuilder report = generateNotAnalysableCPCTReport(0.6, NotAnalysableReason.POST_ANALYSIS_FAIL);
        assertNotNull(report);

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_post_dna_isolation_fail_report.pdf"));
        }
    }

    @NotNull
    private static JasperReportBuilder generateNotAnalysableCPCTReport(final double pathologyTumorEstimate,
            @NotNull final NotAnalysableReason reason) throws IOException {
        final NotAnalysedPatientReport patientReport = ImmutableNotAnalysedPatientReport.of(testSampleReport(pathologyTumorEstimate),
                reason,
                NotAnalysableStudy.CPCT,
                Optional.empty(),
                testBaseReportData().signaturePath());

        return PDFWriter.generateNotAnalysableReport(patientReport);
    }

    @NotNull
    private static SampleReport testSampleReport(final double pathologyTumorPercentage) throws IOException {
        final String sample = "CPCT02991111T";
        return ImmutableSampleReport.of(sample,
                ImmutablePatientTumorLocation.of("CPCT02991111", "Skin", "Melanoma"),
                pathologyTumorPercentage,
                LocalDate.parse("05-Jan-2016", FORMATTER),
                LocalDate.parse("01-Jan-2016", FORMATTER),
                "PREP013V23-QC037V20-SEQ008V25",
                testBaseReportData().centerModel().getAddresseeStringForSample(sample));
    }

    @NotNull
    private static ImmutableEnrichedStructuralVariantLeg.Builder createEnrichedStructuralVariantLegBuilder() {
        return ImmutableEnrichedStructuralVariantLeg.builder().orientation((byte) 1).homology(Strings.EMPTY).position(1);
    }

    @NotNull
    private static ImmutableEnrichedStructuralVariant.Builder createEnrichedStructuralVariantBuilder() {
        return ImmutableEnrichedStructuralVariant.builder().id(Strings.EMPTY).insertSequence(Strings.EMPTY);
    }

    @NotNull
    private static ImmutableEnrichedSomaticVariant.Builder createSomaticVariantBuilder() {
        return SomaticVariantTestBuilderFactory.createEnriched().filter("PASS");
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
}