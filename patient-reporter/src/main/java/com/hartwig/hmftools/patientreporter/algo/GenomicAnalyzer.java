package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.protect.actionability.ClinicalTrialFactory;
import com.hartwig.hmftools.protect.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.protect.bachelor.BachelorData;
import com.hartwig.hmftools.protect.bachelor.BachelorDataLoader;
import com.hartwig.hmftools.protect.chord.ChordDataLoader;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.linx.LinxDataLoader;
import com.hartwig.hmftools.protect.purple.PurpleAnalyzer;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.protect.purple.PurpleDataLoader;
import com.hartwig.hmftools.protect.structural.SvAnalyzer;
import com.hartwig.hmftools.protect.variants.ReportableVariantAnalysis;
import com.hartwig.hmftools.protect.variants.ReportableVariantAnalyzer;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;

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
            @NotNull String chordPredictionTxt) throws IOException {
        PurpleData purpleData =
                PurpleDataLoader.load(tumorSampleId, purpleQCFile, purplePurityTsv, purpleDriverCatalogTsv, purpleSomaticVariantVcf);
        List<EvidenceItem> purpleEvidence = PurpleAnalyzer.run(purpleData, actionabilityAnalyzer, patientPrimaryTumor);

        LinxData linxData = LinxDataLoader.load(linxFusionTsv, linxBreakendTsv, linxViralInsertionTsv, linxDriversTsv);
        List<EvidenceItem> linxEvidence = SvAnalyzer.run(linxData, actionabilityAnalyzer, patientPrimaryTumor);

        BachelorData bachelorData = BachelorDataLoader.load(bachelorTsv, purpleData, linxData, germlineReportingModel);

        ReportableVariantAnalysis reportableVariantAnalysis = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(
                purpleData.somaticVariants(),
                bachelorData.germlineVariants(),
                actionabilityAnalyzer,
                patientPrimaryTumor);

        ChordAnalysis chordAnalysis = ChordDataLoader.load(chordPredictionTxt);

        List<EvidenceItem> allEvidenceItems = Lists.newArrayList();
        allEvidenceItems.addAll(reportableVariantAnalysis.evidenceItems());
        allEvidenceItems.addAll(purpleEvidence);
        allEvidenceItems.addAll(linxEvidence);

        List<EvidenceItem> allEvidenceItemsFiltered = ReportableEvidenceItemFactory.filterBlacklistedEvidence(allEvidenceItems);
        List<EvidenceItem> nonTrials = ReportableEvidenceItemFactory.extractNonTrials(allEvidenceItemsFiltered);

        return ImmutableGenomicAnalysis.builder()
                .impliedPurity(purpleData.purity())
                .hasReliablePurity(purpleData.hasReliablePurity())
                .hasReliableQuality(purpleData.hasReliableQuality())
                .averageTumorPloidy(purpleData.ploidy())
                .tumorSpecificEvidence(nonTrials.stream().filter(EvidenceItem::isOnLabel).collect(Collectors.toList()))
                .clinicalTrials(ClinicalTrialFactory.extractOnLabelTrials(allEvidenceItemsFiltered))
                .offLabelEvidence(nonTrials.stream().filter(item -> !item.isOnLabel()).collect(Collectors.toList()))
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
}
