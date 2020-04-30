package com.hartwig.hmftools.patientreporter;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Locale;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.google.gson.Gson;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.ClinicalTrial;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.EvidenceScope;
import com.hartwig.hmftools.common.actionability.ImmutableClinicalTrial;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ImmutableChordAnalysis;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.ecrf.projections.ImmutablePatientTumorLocation;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.patientreporter.copynumber.ImmutableReportableGainLoss;
import com.hartwig.hmftools.patientreporter.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.patientreporter.homozygousdisruption.ReportableHomozygousDisruption;
import com.hartwig.hmftools.patientreporter.structural.ImmutableReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.variants.ImmutableReportableVariant;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariant;
import com.hartwig.hmftools.patientreporter.viralInsertion.ViralInsertion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class JsonSerializer {

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("dd-MMM-yyyy", Locale.ENGLISH);

    // This code was added by request from Paul R. to generate a COLO829 report in JSON format.
    public static void main(String[] args) {
        Gson gson = new Gson();

        AnalysedPatientReport report = buildCOLO829Report();

        System.out.println(gson.toJson(report));
    }

    @NotNull
    private static AnalysedPatientReport buildCOLO829Report() {
        SampleReport sampleReport = createSkinMelanomaSampleReport("COLO829T");

        String clinicalSummary = "Melanoma sample showing:\n"
                + " - activating BRAF mutation that is associated with response to BRAF-inhibitors (in combination with a MEK-inhibitor)\n"
                + " - complete inactivation of CDKN2A, indicating potential benefit of CDK4/6 inhibitors\n"
                + " - complete inactivation/loss of PTEN likely resulting in an activation of the PI3K-AKT-mTOR pathway "
                + "and indicating potential benefit of mTOR/PI3K inhibitors\n"
                + " - high mutational burden (mutational load (ML) of 180, tumor mutation burden (TMB) of 13.6) that is "
                + "potentially associated with an increased response rate to checkpoint inhibitor immunotherapy";

        List<EvidenceItem> tumorLocationSpecificEvidence = createCOLO829TumorSpecificEvidence();
        List<ClinicalTrial> clinicalTrials = createCOLO829ClinicalTrials();
        List<EvidenceItem> offLabelEvidence = createCOLO829OffLabelEvidence();
        List<ReportableVariant> reportableVariants = createCOLO829SomaticVariants();
        List<ReportableGainLoss> gainsAndLosses = createCOLO829GainsLosses();
        List<ReportableGeneFusion> fusions = Lists.newArrayList();
        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = Lists.newArrayList();
        List<ReportableGeneDisruption> disruptions = createCOLO829Disruptions();
        ChordAnalysis chordAnalysis = createCOLO829ChordAnalysis();
        List<ViralInsertion> viralInsertions = Lists.newArrayList();

        return ImmutableAnalysedPatientReport.builder()
                .sampleReport(sampleReport)
                .impliedPurity(1D)
                .hasReliablePurity(true)
                .hasReliableQuality(true)
                .averageTumorPloidy(3.1)
                .clinicalSummary(clinicalSummary)
                .tumorSpecificEvidence(tumorLocationSpecificEvidence)
                .clinicalTrials(clinicalTrials)
                .offLabelEvidence(offLabelEvidence)
                .reportableVariants(reportableVariants)
                .microsatelliteIndelsPerMb(0.11)
                .tumorMutationalLoad(180)
                .tumorMutationalBurden(13.6)
                .chordAnalysis(chordAnalysis)
                .gainsAndLosses(gainsAndLosses)
                .geneFusions(fusions)
                .geneDisruptions(disruptions)
                .reportableHomozygousDisruptions(reportableHomozygousDisruptions)
                .viralInsertions(viralInsertions)
                .circosPath("path/to/circos.jpg")
                .signaturePath("path/to/signature.jpg")
                .logoRVAPath("path/to/rva_logo.jpg")
                .logoCompanyPath("path/to/company_log.jpg")
                .isCorrectedReport(false)
                .isUnofficialReport(false)
                .comments(Optional.of("This is a test report and is based off COLO829"))
                .build();
    }

    @NotNull
    private static SampleReport createSkinMelanomaSampleReport(@NotNull String sample) {
        SampleMetadata sampleMetadata = ImmutableSampleMetadata.builder()
                .refSampleId(Strings.EMPTY)
                .refSampleBarcode("FR12345678")
                .tumorSampleId(sample)
                .tumorSampleBarcode("FR23456789")
                .build();

        return ImmutableSampleReport.builder()
                .sampleMetadata(sampleMetadata)
                .patientTumorLocation(ImmutablePatientTumorLocation.of(Strings.EMPTY, "Skin", "Melanoma"))
                .refArrivalDate(LocalDate.parse("01-Jan-2020", DATE_FORMATTER))
                .tumorArrivalDate(LocalDate.parse("05-Jan-2020", DATE_FORMATTER))
                .purityShallowSeq(Strings.EMPTY)
                .labProcedures("PREP013V23-QC037V20-SEQ008V25")
                .requesterName("Paul")
                .requesterEmail("paul@hartwig.com")
                .addressee("HMF Testing Center")
                .hospitalName(Strings.EMPTY)
                .hospitalPIName(Strings.EMPTY)
                .hospitalPIEmail(Strings.EMPTY)
                .cohort("TEST")
                .projectName("TEST")
                .submissionId("10")
                .hospitalPatientId("4567")
                .hospitalPathologySampleId("1234")
                .build();
    }

    @NotNull
    private static List<EvidenceItem> createCOLO829TumorSpecificEvidence() {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();

        ImmutableEvidenceItem.Builder onLabelBuilder = evidenceBuilder().isOnLabel(true);

        evidenceItems.add(onLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Binimetinib + Encorafenib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("V600E")
                .source(ActionabilitySource.ONCOKB)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(onLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Cobimetinib + Vemurafenib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("V600E")
                .source(ActionabilitySource.ONCOKB)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(onLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Dabrafenib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("V600E")
                .source(ActionabilitySource.ONCOKB)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(onLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Dabrafenib + Trametinib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("V600E")
                .source(ActionabilitySource.ONCOKB)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(onLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Trametinib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("V600E")
                .source(ActionabilitySource.ONCOKB)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(onLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Vemurafenib")
                .level(EvidenceLevel.LEVEL_A)
                .response("Responsive")
                .reference("V600E")
                .source(ActionabilitySource.ONCOKB)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(onLabelBuilder.event("BRAF p.Val600Glu")
                .drug("RO4987655")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:208")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.GENE_LEVEL)
                .build());

        return evidenceItems;
    }

    @NotNull
    private static List<ClinicalTrial> createCOLO829ClinicalTrials() {
        List<ClinicalTrial> trials = Lists.newArrayList();
        ImmutableClinicalTrial.Builder iclusionBuilder =
                ImmutableClinicalTrial.builder().cancerType(Strings.EMPTY).isOnLabel(true).source(ActionabilitySource.ICLUSION);

        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.GENE_LEVEL)
                .acronym("DRUP")
                .reference("EXT10299 (NL54757.031.16)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.GENE_LEVEL)
                .acronym("EBIN (EORTC-1612-MG)")
                .reference("EXT11284 (NL67202.031.18)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.GENE_LEVEL)
                .acronym("LXH254 in tumors with MAPK pathway alterations")
                .reference("EXT10453 (NL55506.078.15)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.GENE_LEVEL)
                .acronym("POLARIS")
                .reference("EXT11388 (NL69569.028.19)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.SPECIFIC)
                .acronym("PROCLAIM-001")
                .reference("EXT10241 (NL59299.042.17)")
                .build());
        trials.add(iclusionBuilder.event("BRAF p.Val600Glu")
                .scope(EvidenceScope.SPECIFIC)
                .acronym("REDUCTOR")
                .reference("EXT6690 (NL45261.031.13)")
                .build());
        trials.add(iclusionBuilder.event("CDKN2A p.Ala68fs")
                .scope(EvidenceScope.GENE_LEVEL)
                .acronym("DRUP")
                .reference("EXT10299 (NL54757.031.16)")
                .build());

        return trials;
    }

    @NotNull
    private static List<EvidenceItem> createCOLO829OffLabelEvidence() {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();

        ImmutableEvidenceItem.Builder offLabelBuilder = evidenceBuilder().isOnLabel(false);

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Alpelisib + Cetuximab + Encorafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.GENE_LEVEL)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Bevacizumab")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("CI-1040")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Cetuximab")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("BRAF:V600E")
                .source(ActionabilitySource.CGI)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Cetuximab + Encorafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.GENE_LEVEL)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Cetuximab + Irinotecan + Vemurafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Dabrafenib + Panitumumab + Trametinib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Irinotecan")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Oxaliplatin")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("variant:12")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Panitumumab")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("BRAF:V600E")
                .source(ActionabilitySource.CGI)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("BRAF p.Val600Glu")
                .drug("Vemurafenib")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("variant:17")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.GENE_LEVEL)
                .build());

        evidenceItems.add(offLabelBuilder.event("PTEN Deletion")
                .drug("EGFR mAB inhibitor")
                .level(EvidenceLevel.LEVEL_B)
                .response("Resistant")
                .reference("PTEN:del")
                .source(ActionabilitySource.CGI)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        evidenceItems.add(offLabelBuilder.event("PTEN Deletion")
                .drug("Everolimus")
                .level(EvidenceLevel.LEVEL_B)
                .response("Responsive")
                .reference("variant:213")
                .source(ActionabilitySource.CIVIC)
                .scope(EvidenceScope.SPECIFIC)
                .build());

        return evidenceItems;
    }

    @NotNull
    private static List<ReportableVariant> createCOLO829SomaticVariants() {
        ReportableVariant variant1 = ImmutableReportableVariant.builder()
                .gene("BRAF")
                .position(140453136)
                .chromosome("7")
                .ref("T")
                .alt("A")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.1799T>A")
                .canonicalHgvsProteinImpact("p.Val600Glu")
                .notifyClinicalGeneticist(false)
                .driverCategory(DriverCategory.ONCO)
                .gDNA("7:140453136")
                .alleleReadCount(154)
                .totalReadCount(225)
                .allelePloidy(4.1)
                .totalPloidy(6.1)
                .hotspot(Hotspot.HOTSPOT)
                .biallelic(false)
                .driverLikelihood(1D)
                .clonalLikelihood(1D)
                .build();

        ReportableVariant variant2 = ImmutableReportableVariant.builder()
                .gene("CDKN2A")
                .position(21971153)
                .chromosome("9")
                .ref("CCG")
                .alt("C")
                .canonicalCodingEffect(CodingEffect.NONSENSE_OR_FRAMESHIFT)
                .canonicalHgvsCodingImpact("c.203_204delCG")
                .canonicalHgvsProteinImpact("p.Ala68fs")
                .notifyClinicalGeneticist(false)
                .driverCategory(DriverCategory.TSG)
                .gDNA("9:21971153")
                .alleleReadCount(95)
                .totalReadCount(95)
                .allelePloidy(2D)
                .totalPloidy(2D)
                .hotspot(Hotspot.NEAR_HOTSPOT)
                .biallelic(true)
                .clonalLikelihood(1D)
                .driverLikelihood(0.9)
                .build();

        ReportableVariant variant3 = ImmutableReportableVariant.builder()
                .gene("TERT")
                .position(1295228)
                .chromosome("5")
                .ref("GG")
                .alt("AA")
                .canonicalCodingEffect(CodingEffect.NONE)
                .canonicalHgvsCodingImpact("c.-125_-124delCCinsTT")
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .notifyClinicalGeneticist(false)
                .driverCategory(DriverCategory.ONCO)
                .gDNA("5:1295228")
                .alleleReadCount(49)
                .totalReadCount(49)
                .allelePloidy(2D)
                .totalPloidy(2D)
                .hotspot(Hotspot.HOTSPOT)
                .biallelic(false)
                .clonalLikelihood(1D)
                .driverLikelihood(0.85)
                .build();

        ReportableVariant variant4 = ImmutableReportableVariant.builder()
                .gene("SF3B1")
                .position(198266779)
                .chromosome("2")
                .ref("C")
                .alt("T")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.2153C>T")
                .canonicalHgvsProteinImpact("p.Pro718Leu")
                .notifyClinicalGeneticist(false)
                .driverCategory(DriverCategory.ONCO)
                .gDNA("2:198266779")
                .alleleReadCount(76)
                .totalReadCount(115)
                .allelePloidy(2D)
                .totalPloidy(3.1)
                .hotspot(Hotspot.NON_HOTSPOT)
                .biallelic(false)
                .clonalLikelihood(1D)
                .driverLikelihood(0.5)
                .build();

        ReportableVariant variant5 = ImmutableReportableVariant.builder()
                .gene("TP63")
                .position(189604330)
                .chromosome("3")
                .ref("G")
                .alt("T")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("c.1497G>T")
                .canonicalHgvsProteinImpact("p.Met499Ile")
                .notifyClinicalGeneticist(false)
                .driverCategory(DriverCategory.TSG)
                .gDNA("3:189604330")
                .alleleReadCount(52)
                .totalReadCount(119)
                .allelePloidy(1.8)
                .totalPloidy(4D)
                .hotspot(Hotspot.NON_HOTSPOT)
                .biallelic(false)
                .clonalLikelihood(1D)
                .driverLikelihood(0.1)
                .build();

        return Lists.newArrayList(variant1, variant2, variant3, variant4, variant5);
    }

    @NotNull
    private static List<ReportableGainLoss> createCOLO829GainsLosses() {
        ReportableGainLoss gainLoss1 = ImmutableReportableGainLoss.builder()
                .chromosome("10")
                .chromosomeBand("q23.31")
                .gene("PTEN")
                .copies(0)
                .interpretation(CopyNumberInterpretation.PARTIAL_LOSS)
                .build();

        return Lists.newArrayList(gainLoss1);
    }

    @NotNull
    private static ChordAnalysis createCOLO829ChordAnalysis() {
        double brca1Value = 0D;
        double brca2Value = 0D;

        return ImmutableChordAnalysis.builder()
                .noneValue(1 - (brca1Value + brca2Value))
                .BRCA1Value(brca1Value)
                .BRCA2Value(brca2Value)
                .hrdValue(brca1Value + brca2Value)
                .predictedResponseValue(brca1Value + brca2Value > 0.5)
                .build();
    }

    @NotNull
    private static List<ReportableGeneDisruption> createCOLO829Disruptions() {
        ReportableGeneDisruption disruption1 = createDisruptionBuilder().location("10q23.31")
                .gene("PTEN")
                .range("Intron 5 -> Intron 6")
                .type("DEL")
                .ploidy(2D)
                .undisruptedCopyNumber(0)
                .build();

        return Lists.newArrayList(disruption1);
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder evidenceBuilder() {
        return ImmutableEvidenceItem.builder().drugsType(Strings.EMPTY).cancerType(Strings.EMPTY);
    }

    @NotNull
    private static ImmutableReportableGeneDisruption.Builder createDisruptionBuilder() {
        return ImmutableReportableGeneDisruption.builder().firstAffectedExon(1);
    }
}
