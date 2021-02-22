package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceFile;
import com.hartwig.hmftools.patientreporter.actionability.ClinicalTrialFactory;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.protect.bachelor.BachelorData;
import com.hartwig.hmftools.protect.bachelor.BachelorDataLoader;
import com.hartwig.hmftools.protect.chord.ChordDataLoader;
import com.hartwig.hmftools.protect.germline.GermlineReportingModel;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.linx.LinxDataLoader;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.protect.purple.PurpleDataLoader;
import com.hartwig.hmftools.protect.purple.ReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariantFactory;

import org.jetbrains.annotations.NotNull;

public class GenomicAnalyzer {

    @NotNull
    private final GermlineReportingModel germlineReportingModel;

    public GenomicAnalyzer(@NotNull final GermlineReportingModel germlineReportingModel) {
        this.germlineReportingModel = germlineReportingModel;
    }

    @NotNull
    public GenomicAnalysis run(@NotNull String tumorSampleId, @NotNull String purplePurityTsv, @NotNull String purpleQCFile,
            @NotNull String purpleDriverCatalogSomaticTsv, @NotNull String purpleDriverCatalogGermlineTsv,
            @NotNull String purpleSomaticVariantVcf, @NotNull String purpleGermlineVariantVcf, @NotNull String bachelorTsv, @NotNull String linxFusionTsv,
            @NotNull String linxBreakendTsv, @NotNull String linxViralInsertionTsv, @NotNull String linxDriversTsv,
            @NotNull String chordPredictionTxt, @NotNull String protectEvidenceTsv) throws IOException {
        PurpleData purpleData = PurpleDataLoader.load(tumorSampleId,
                purpleQCFile,
                purplePurityTsv,
                purpleDriverCatalogSomaticTsv,
                purpleSomaticVariantVcf,
                purpleDriverCatalogGermlineTsv,
                purpleGermlineVariantVcf);

        LinxData linxData = LinxDataLoader.load(linxFusionTsv, linxBreakendTsv, linxViralInsertionTsv, linxDriversTsv);

        List<DriverCatalog> germlineDriverCatalog = DriverCatalogFile.read(purpleDriverCatalogGermlineTsv);

        BachelorData bachelorData = BachelorDataLoader.load(bachelorTsv, purpleData, linxData, germlineReportingModel);

        ReportableVariantAnalysis reportableVariantAnalysis =
                mergeSomaticAndGermlineVariants(purpleData.somaticVariants(), bachelorData.germlineVariants());

        ChordAnalysis chordAnalysis = ChordDataLoader.load(chordPredictionTxt);

        List<ProtectEvidence> reportableEvidenceItems = extractReportableEvidenceItems(protectEvidenceTsv);

        List<ProtectEvidence> nonTrials = ReportableEvidenceItemFactory.extractNonTrials(reportableEvidenceItems);

        return ImmutableGenomicAnalysis.builder()
                .impliedPurity(purpleData.purity())
                .hasReliablePurity(purpleData.hasReliablePurity())
                .hasReliableQuality(purpleData.hasReliableQuality())
                .averageTumorPloidy(purpleData.ploidy())
                .tumorSpecificEvidence(nonTrials.stream().filter(ProtectEvidence::onLabel).collect(Collectors.toList()))
                .clinicalTrials(ClinicalTrialFactory.extractOnLabelTrials(reportableEvidenceItems))
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
    private static ReportableVariantAnalysis mergeSomaticAndGermlineVariants(@NotNull List<ReportableVariant> reportableSomaticVariants,
            @NotNull List<ReportableVariant> reportableGermlineVariants) {
        List<ReportableVariant> allReportableVariants =
                ReportableVariantFactory.mergeVariantLists(reportableGermlineVariants, reportableSomaticVariants);

        return ImmutableReportableVariantAnalysis.builder().variantsToReport(allReportableVariants).build();
    }

    @NotNull
    private static List<ProtectEvidence> extractReportableEvidenceItems(@NotNull String protectEvidenceTsv) throws IOException {
        List<ProtectEvidence> evidences = ProtectEvidenceFile.read(protectEvidenceTsv);

        List<ProtectEvidence> reportableEvidenceItems = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (evidence.reported()) {
                reportableEvidenceItems.add(evidence);
            }
        }
        return reportableEvidenceItems;
    }
}
