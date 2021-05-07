package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceFile;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakendQCStatus;
import com.hartwig.hmftools.patientreporter.actionability.ClinicalTrialFactory;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakendFactory;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakend;
import com.hartwig.hmftools.patientreporter.germline.GermlineCondition;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingEntry;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusDbModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusSummaryModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.ImmutableReportableVirusBreakendTotal;
import com.hartwig.hmftools.patientreporter.virusbreakend.ImmutableReportableVirusbreakend;
import com.hartwig.hmftools.patientreporter.virusbreakend.ReportableVirusBreakendTotal;
import com.hartwig.hmftools.patientreporter.virusbreakend.ReportableVirusbreakend;
import com.hartwig.hmftools.protect.chord.ChordDataLoader;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.linx.LinxDataLoader;
import com.hartwig.hmftools.protect.purple.ReportableVariantSource;
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
            @NotNull String linxBreakendTsv, @NotNull String linxDriversTsv, @NotNull String chordPredictionTxt,
            @NotNull String protectEvidenceTsv, @NotNull String virusBreakendTsv, @NotNull String peachgenotypeTsv,
            @NotNull GermlineReportingModel germlineReportingModel, @NotNull LimsGermlineReportingLevel germlineReportingLevel,
            @NotNull VirusDbModel virusDbModel, @NotNull VirusSummaryModel virusSummaryModel) throws IOException {
        PurpleData purpleData = PurpleDataLoader.load(tumorSampleId,
                purpleQCFile,
                purplePurityTsv,
                purpleDriverCatalogSomaticTsv,
                purpleSomaticVariantVcf,
                purpleDriverCatalogGermlineTsv,
                purpleGermlineVariantVcf);

        LinxData linxData = LinxDataLoader.load(linxFusionTsv, linxBreakendTsv, linxDriversTsv);

        List<VirusBreakend> virusBreakends = VirusBreakendFactory.readVirusBreakend(virusBreakendTsv);

        List<ReportableVirusbreakend> virusBreakendsReportable = Lists.newArrayList();
        Set<String> postiveSummary = Sets.newHashSet();
        Set<String> negativeSummary = Sets.newHashSet();

        Set<String> summary = Sets.newHashSet();
        for (VirusBreakend virusBreakend : virusBreakends) {
            if (virusBreakend.QCStatus() != VirusBreakendQCStatus.LOW_VIRAL_COVERAGE) {
                if (virusBreakend.integrations() >= 1) {

                    String virusName = virusDbModel.findVirus(virusBreakend.referenceTaxid());
                    virusBreakendsReportable.add(ImmutableReportableVirusbreakend.builder()
                            .virusName(virusName)
                            .integrations(virusBreakend.integrations())
                            .build());

                    if (virusSummaryModel.mapIdtoVirusName(virusBreakend.taxidSpecies())) {
                        postiveSummary.add(virusSummaryModel.findVirusSummary(virusBreakend.taxidSpecies()) + " positive");
                    }
                }
            }
        }

        for (VirusBreakend virusBreakend : virusBreakends) {
            for (String virus : virusSummaryModel.virussen()) {
                if (!postiveSummary.contains(virus)) {
                    negativeSummary.add(virusSummaryModel.findVirusSummary(virusBreakend.taxidSpecies()) + " negative");
                }
            }
        }

        summary.addAll(postiveSummary);
        summary.addAll(negativeSummary);

        ReportableVirusBreakendTotal reportableVirusBreakendTotal = ImmutableReportableVirusBreakendTotal.builder()
                .reportableVirussen(virusBreakendsReportable)
                .virusNameSummary(summary.toString())
                .build();

        List<PeachGenotype> peachGenotypes = PeachGenotypeFile.read(peachgenotypeTsv);

        List<ReportableVariant> reportableVariants =
                ReportableVariantFactory.mergeVariantLists(purpleData.germlineVariants(), purpleData.somaticVariants());

        List<ReportableVariantNotify> reportableVariantsWithNotify =
                determineNotify(reportableVariants, germlineReportingModel, germlineReportingLevel);

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
                .virusBreakends(reportableVirusBreakendTotal)
                .peachGenotypes(peachGenotypes)
                .build();
    }

    @NotNull
    private static List<ReportableVariantNotify> determineNotify(List<ReportableVariant> reportableVariants,
            @NotNull GermlineReportingModel germlineReportingModel, @NotNull LimsGermlineReportingLevel germlineReportingLevel) {
        List<ReportableVariantNotify> reportableVariantNotifies = Lists.newArrayList();

        boolean notify;
        for (ReportableVariant reportableVariant : reportableVariants) {
            if (reportableVariant.source() == ReportableVariantSource.GERMLINE) {
                notify = notifyAboutVariant(reportableVariant, germlineReportingModel, germlineReportingLevel);

            } else {
                notify = false;

            }
            reportableVariantNotifies.add(ImmutableReportableVariantNotify.builder()
                    .reportableVariant(reportableVariant)
                    .notifyVariant(notify)
                    .build());
        }

        return reportableVariantNotifies;
    }

    private static boolean notifyAboutVariant(@NotNull ReportableVariant variant, @NotNull GermlineReportingModel germlineReportingModel,
            @NotNull LimsGermlineReportingLevel germlineReportingLevel) {
        boolean notifyVariant = false;
        if (variant.source() == ReportableVariantSource.GERMLINE) {
            GermlineReportingEntry reportingEntry = germlineReportingModel.entryForGene(variant.gene());
            if (reportingEntry != null) {
                if (reportingEntry.notifyClinicalGeneticist() == GermlineCondition.ONLY_GERMLINE_HOM) {
                    String conditionFilter = reportingEntry.conditionFilter();
                    if (conditionFilter != null) {
                        notifyVariant = variant.genotypeStatus().simplifiedDisplay().equals(conditionFilter);
                    }

                } else if (reportingEntry.notifyClinicalGeneticist() == GermlineCondition.ONLY_SPECIFIC_VARIANT) {
                    String conditionFilter = reportingEntry.conditionFilter();
                    if (conditionFilter != null) {
                        notifyVariant = variant.canonicalHgvsProteinImpact().equals(conditionFilter);
                    }
                } else if (reportingEntry.notifyClinicalGeneticist() == GermlineCondition.ALWAYS) {
                    notifyVariant = true;
                }
            }
        }
        return notifyVariant && germlineReportingModel.notifyAboutGene(variant.gene(), germlineReportingLevel);
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
