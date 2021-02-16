package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceFile;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.patientreporter.actionability.ClinicalTrialFactory;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.protect.bachelor.BachelorData;
import com.hartwig.hmftools.protect.bachelor.BachelorDataLoader;
import com.hartwig.hmftools.protect.chord.ChordDataLoader;
import com.hartwig.hmftools.protect.germline.GermlineReportingModel;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.linx.LinxDataLoader;
import com.hartwig.hmftools.protect.purple.DriverInterpretation;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.protect.purple.PurpleDataLoader;
import com.hartwig.hmftools.protect.purple.ReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariantFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GenomicAnalyzer {

    @NotNull
    private final ActionabilityAnalyzer actionabilityAnalyzer;
    @NotNull
    private final GermlineReportingModel germlineReportingModel;

    public GenomicAnalyzer(@NotNull final ActionabilityAnalyzer actionabilityAnalyzer,
            @NotNull final GermlineReportingModel germlineReportingModel) {
        this.actionabilityAnalyzer = actionabilityAnalyzer;
        this.germlineReportingModel = germlineReportingModel;
    }

    @NotNull
    public GenomicAnalysis run(@NotNull String tumorSampleId, @Nullable PatientPrimaryTumor patientPrimaryTumor,
            @NotNull String purplePurityTsv, @NotNull String purpleQCFile, @NotNull String purpleDriverCatalogTsv,
            @NotNull String purpleSomaticVariantVcf, @NotNull String bachelorTsv, @NotNull String linxFusionTsv,
            @NotNull String linxBreakendTsv, @NotNull String linxViralInsertionTsv, @NotNull String linxDriversTsv,
            @NotNull String chordPredictionTxt, @NotNull String protectEvidenceTsv) throws IOException {
        PurpleData purpleData =
                PurpleDataLoader.load(tumorSampleId, purpleQCFile, purplePurityTsv, purpleDriverCatalogTsv, purpleSomaticVariantVcf);
        List<EvidenceItem> purpleEvidence = determinePurpleEvidence(purpleData, actionabilityAnalyzer, patientPrimaryTumor);

        LinxData linxData = LinxDataLoader.load(linxFusionTsv, linxBreakendTsv, linxViralInsertionTsv, linxDriversTsv);
        List<EvidenceItem> linxEvidence = determineLinxEvidence(linxData, actionabilityAnalyzer, patientPrimaryTumor);

        BachelorData bachelorData = BachelorDataLoader.load(bachelorTsv, purpleData, linxData, germlineReportingModel);

        ReportableVariantAnalysis reportableVariantAnalysis = mergeSomaticAndGermlineVariants(
                purpleData.somaticVariants(),
                bachelorData.germlineVariants(),
                actionabilityAnalyzer,
                patientPrimaryTumor);

        ChordAnalysis chordAnalysis = ChordDataLoader.load(chordPredictionTxt);

        List<ProtectEvidence> allEvidenceItems = Lists.newArrayList();

        List<ProtectEvidence> evidences = ProtectEvidenceFile.read(protectEvidenceTsv);
        for (ProtectEvidence evidence: evidences) {
            if (evidence.reported()) {
                allEvidenceItems.add(evidence);
            }
        }



      //  List<EvidenceItem> allEvidenceItemsFiltered = ReportableEvidenceItemFactory.filterBlacklistedEvidence(allEvidenceItems);
        List<ProtectEvidence> nonTrials = ReportableEvidenceItemFactory.extractNonTrials(allEvidenceItems);

        return ImmutableGenomicAnalysis.builder()
                .impliedPurity(purpleData.purity())
                .hasReliablePurity(purpleData.hasReliablePurity())
                .hasReliableQuality(purpleData.hasReliableQuality())
                .averageTumorPloidy(purpleData.ploidy())
                .tumorSpecificEvidence(nonTrials.stream().filter(ProtectEvidence::onLabel).collect(Collectors.toList()))
                .clinicalTrials(ClinicalTrialFactory.extractOnLabelTrials(allEvidenceItems))
                .offLabelEvidence(nonTrials.stream().filter(item -> !item.onLabel()).collect(Collectors.toList()))
                .reportableVariants(reportableVariantAnalysis.variantsToReport())
                .microsatelliteIndelsPerMb(purpleData.microsatelliteIndelsPerMb())
                .microsatelliteStatus(purpleData.microsatelliteStatus())
                .tumorMutationalLoad(purpleData.tumorMutationalLoad())
                .tumorMutationalLoadStatus(purpleData.tumorMutationalLoadStatus())
                .tumorMutationalBurden(purpleData.tumorMutationalBurdenPerMb())
                .chordHrdValue(chordAnalysis.hrdValue())
                .chordHrdStatus(chordAnalysis.hrStatus())
                .gainsAndLosses(purpleData.copyNumberAlterations())
                .geneFusions(linxData.fusions())
                .geneDisruptions(linxData.geneDisruptions())
                .homozygousDisruptions(linxData.homozygousDisruptions())
                .viralInsertions(linxData.viralInsertions())
                .build();
    }

    @NotNull
    private static List<EvidenceItem> determinePurpleEvidence(@NotNull PurpleData purpleData,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer, @Nullable PatientPrimaryTumor patientPrimaryTumor) {
        String primaryTumorLocation = patientPrimaryTumor != null ? patientPrimaryTumor.location() : null;
        Map<ReportableGainLoss, List<EvidenceItem>> evidencePerGeneCopyNumber =
                actionabilityAnalyzer.evidenceForCopyNumbers(purpleData.copyNumberAlterations(), primaryTumorLocation);

        return ReportableEvidenceItemFactory.toReportableFlatList(evidencePerGeneCopyNumber);
    }

    @NotNull
    private static List<EvidenceItem> determineLinxEvidence(@NotNull LinxData linxData, @NotNull ActionabilityAnalyzer actionabilityAnalyzer,
            @Nullable PatientPrimaryTumor patientPrimaryTumor) {
        String primaryTumorLocation = patientPrimaryTumor != null ? patientPrimaryTumor.location() : null;
        Map<LinxFusion, List<EvidenceItem>> evidencePerFusion =
                actionabilityAnalyzer.evidenceForFusions(linxData.fusions(), primaryTumorLocation);

        return ReportableEvidenceItemFactory.toReportableFlatList(evidencePerFusion);
    }

    @NotNull
    private static ReportableVariantAnalysis mergeSomaticAndGermlineVariants(@NotNull List<ReportableVariant> reportableSomaticVariants,
            @NotNull List<ReportableVariant> reportableGermlineVariants, @NotNull ActionabilityAnalyzer actionabilityAnalyzer,
            @Nullable PatientPrimaryTumor patientPrimaryTumor) {
        List<ReportableVariant> allReportableVariants =
                ReportableVariantFactory.mergeVariantLists(reportableGermlineVariants, reportableSomaticVariants);

        String primaryTumorLocation = patientPrimaryTumor != null ? patientPrimaryTumor.location() : null;
        // Extract somatic evidence for high drivers variants only (See DEV-824)
        Map<ReportableVariant, List<EvidenceItem>> evidencePerVariant =
                filterHighDriverLikelihood(actionabilityAnalyzer.evidenceForAllVariants(allReportableVariants, primaryTumorLocation));

        return ImmutableReportableVariantAnalysis.builder()
                .variantsToReport(allReportableVariants)
                .evidenceItems(ReportableEvidenceItemFactory.toReportableFlatList(evidencePerVariant))
                .build();
    }

    @NotNull
    private static Map<ReportableVariant, List<EvidenceItem>> filterHighDriverLikelihood(
            @NotNull Map<? extends Variant, List<EvidenceItem>> evidenceForAllVariants) {
        Map<ReportableVariant, List<EvidenceItem>> evidencePerHighDriverVariant = Maps.newHashMap();
        for (Map.Entry<? extends Variant, List<EvidenceItem>> entry : evidenceForAllVariants.entrySet()) {
            ReportableVariant variant = (ReportableVariant) entry.getKey();
            if (DriverInterpretation.interpret(variant.driverLikelihood()) == DriverInterpretation.HIGH) {
                evidencePerHighDriverVariant.put(variant, entry.getValue());
            }
        }
        return evidencePerHighDriverVariant;
    }
}
