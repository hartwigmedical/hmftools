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
import com.hartwig.hmftools.common.virusbreakend.VirusBreakend;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakendFactory;
import com.hartwig.hmftools.patientreporter.PatientReporterConfig;
import com.hartwig.hmftools.patientreporter.actionability.ClinicalTrialFactory;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.germline.GermlineCondition;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingEntry;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.ReportableVirusBreakendTotal;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusBlacklistModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusBreakendReportableFactory;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusDbModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusSummaryModel;
import com.hartwig.hmftools.protect.chord.ChordDataLoader;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.linx.LinxDataLoader;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.protect.purple.PurpleDataLoader;
import com.hartwig.hmftools.protect.purple.ReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariantFactory;
import com.hartwig.hmftools.protect.purple.ReportableVariantSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GenomicAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(GenomicAnalyzer.class);

    @NotNull
    private final GermlineReportingModel germlineReportingModel;
    @NotNull
    private final VirusDbModel virusDbModel;
    @NotNull
    private final VirusSummaryModel virusSummaryModel;
    @NotNull
    private final VirusBlacklistModel virusBlackListModel;

    public GenomicAnalyzer(@NotNull final GermlineReportingModel germlineReportingModel, @NotNull final VirusDbModel virusDbModel,
            @NotNull final VirusSummaryModel virusSummaryModel, @NotNull final VirusBlacklistModel virusBlackListModel) {
        this.germlineReportingModel = germlineReportingModel;
        this.virusDbModel = virusDbModel;
        this.virusSummaryModel = virusSummaryModel;
        this.virusBlackListModel = virusBlackListModel;
    }

    @NotNull
    public GenomicAnalysis run(@NotNull String tumorSampleId, @NotNull PatientReporterConfig config,
            @NotNull LimsGermlineReportingLevel germlineReportingLevel) throws IOException {
        PurpleData purpleData = PurpleDataLoader.load(tumorSampleId,
                config.purpleQcFile(),
                config.purplePurityTsv(),
                config.purpleSomaticDriverCatalogTsv(),
                config.purpleSomaticVariantVcf(),
                config.purpleGermlineDriverCatalogTsv(),
                config.purpleGermlineVariantVcf(),
                config.purpleSomaticCopyNumberTsv());

        LinxData linxData = LinxDataLoader.load(config.linxFusionTsv(), config.linxBreakendTsv(), config.linxDriverCatalogTsv());

        List<VirusBreakend> virusBreakends = VirusBreakendFactory.readVirusBreakend(config.virusBreakendTsv());

        ReportableVirusBreakendTotal reportableVirusBreakendTotal =
                VirusBreakendReportableFactory.analyzeVirusBreakend(virusBreakends, virusDbModel, virusSummaryModel, virusBlackListModel);

        List<PeachGenotype> peachGenotypes = PeachGenotypeFile.read(config.peachGenotypeTsv());

        List<ReportableVariant> reportableVariants =
                ReportableVariantFactory.mergeVariantLists(purpleData.germlineVariants(), purpleData.somaticVariants());

        List<ReportableVariantNotify> reportableVariantsWithNotify =
                determineNotify(reportableVariants, germlineReportingModel, germlineReportingLevel);

        ChordAnalysis chordAnalysis = ChordDataLoader.load(config.chordPredictionTxt());

        List<ProtectEvidence> reportableEvidenceItems = extractReportableEvidenceItems(config.protectEvidenceTsv());

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
                .cnPerChromosome(purpleData.cnPerChromosome())
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
        Set<String> germlineGenes = Sets.newHashSet();
        for (ReportableVariant reportableVariant : reportableVariants) {
            if (reportableVariant.source() == ReportableVariantSource.GERMLINE && reportableVariant.localPhaseSet() == null) {
                germlineGenes.add(reportableVariant.gene());
            }
        }

        boolean notify;
        for (ReportableVariant reportableVariant : reportableVariants) {
            if (reportableVariant.source() == ReportableVariantSource.GERMLINE) {
                notify = notifyAboutVariant(reportableVariant, germlineReportingModel, germlineReportingLevel, germlineGenes);

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
            @NotNull LimsGermlineReportingLevel germlineReportingLevel, @NotNull Set<String> germlineGenes) {
        boolean notifyVariant = false;
        if (variant.source() == ReportableVariantSource.GERMLINE) {
            GermlineReportingEntry reportingEntry = germlineReportingModel.entryForGene(variant.gene());
            if (reportingEntry != null) {
                if (reportingEntry.notifyClinicalGeneticist() == GermlineCondition.ONLY_GERMLINE_HOM) {
                    String conditionFilter = reportingEntry.conditionFilter();
                    if (conditionFilter != null) {
                        if (germlineGenes.contains(variant.gene())) {
                            notifyVariant = true;
                        } else {
                            notifyVariant = variant.genotypeStatus().simplifiedDisplay().equals(conditionFilter);
                        }
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
