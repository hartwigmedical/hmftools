package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceFile;
import com.hartwig.hmftools.patientreporter.actionability.ClinicalTrialFactory;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.viralbreakend.ViralBreakendFactory;
import com.hartwig.hmftools.patientreporter.viralbreakend.Viralbreakend;
import com.hartwig.hmftools.protect.chord.ChordDataLoader;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.linx.LinxDataLoader;
import com.hartwig.hmftools.protect.linx.ViralInsertionAnalyzer;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.protect.purple.PurpleDataLoader;
import com.hartwig.hmftools.protect.purple.ReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GenomicAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(GenomicAnalyzer.class);

    public GenomicAnalyzer() {
    }

    @NotNull
    public GenomicAnalysis run(@NotNull String tumorSampleId, @NotNull String purplePurityTsv, @NotNull String purpleQCFile,
            @NotNull String purpleDriverCatalogSomaticTsv, @NotNull String purpleDriverCatalogGermlineTsv,
            @NotNull String purpleSomaticVariantVcf, @NotNull String purpleGermlineVariantVcf, @NotNull String linxFusionTsv,
            @NotNull String linxBreakendTsv, @NotNull String linxDriversTsv,
            @NotNull String chordPredictionTxt, @NotNull String protectEvidenceTsv, @NotNull String viralBreakendTsv) throws IOException {
        PurpleData purpleData = PurpleDataLoader.load(tumorSampleId,
                purpleQCFile,
                purplePurityTsv,
                purpleDriverCatalogSomaticTsv,
                purpleSomaticVariantVcf,
                purpleDriverCatalogGermlineTsv,
                purpleGermlineVariantVcf);

        LinxData linxData = LinxDataLoader.load(linxFusionTsv, linxBreakendTsv, linxDriversTsv);

        List<Viralbreakend> viralBreakends = ViralBreakendFactory.readViralBreakend(viralBreakendTsv);

        List<ReportableVariant> reportableVariants =
                ReportableVariantFactory.mergeVariantLists(purpleData.germlineVariants(), purpleData.somaticVariants());

        ChordAnalysis chordAnalysis = ChordDataLoader.load(chordPredictionTxt);

        List<ProtectEvidence> reportableEvidenceItems = extractReportableEvidenceItems(protectEvidenceTsv);

        List<ProtectEvidence> nonTrialsOnLabel = ReportableEvidenceItemFactory.extractNonTrialsOnLabel(reportableEvidenceItems);
        List<ProtectEvidence> clinicalTrialsOnLabel = ClinicalTrialFactory.extractOnLabelTrials(reportableEvidenceItems);
        List<ProtectEvidence> nonTrialsOffLabel = ReportableEvidenceItemFactory.extractNonTrialsOffLabel(reportableEvidenceItems);

        return ImmutableGenomicAnalysis.builder()
                .impliedPurity(purpleData.purity())
                .hasReliablePurity(purpleData.hasReliablePurity())
                .hasReliableQuality(purpleData.hasReliableQuality())
                .averageTumorPloidy(purpleData.ploidy())
                .tumorSpecificEvidence(nonTrialsOnLabel)
                .clinicalTrials(clinicalTrialsOnLabel)
                .offLabelEvidence(nonTrialsOffLabel)
                .reportableVariants(reportableVariants)
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
                .viralBreakends(viralBreakends)
                .build();
    }

    @NotNull
    private static List<ProtectEvidence> extractReportableEvidenceItems(@NotNull String protectEvidenceTsv) throws IOException {
        LOGGER.info("Loading PROTECT data from {}", protectEvidenceTsv);
        List<ProtectEvidence> evidences = ProtectEvidenceFile.read(protectEvidenceTsv);

        List<ProtectEvidence> reportableEvidenceItems = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (evidence.reported()) {
                reportableEvidenceItems.add(evidence);
            }
        }
        LOGGER.info(" Loaded {} reportable evidence items", reportableEvidenceItems.size());
        return reportableEvidenceItems;
    }
}
