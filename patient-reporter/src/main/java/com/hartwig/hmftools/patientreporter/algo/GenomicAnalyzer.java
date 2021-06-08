package com.hartwig.hmftools.patientreporter.algo;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceFile;
import com.hartwig.hmftools.patientreporter.PatientReporterConfig;
import com.hartwig.hmftools.patientreporter.actionability.ClinicalTrialFactory;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingModel;
import com.hartwig.hmftools.protect.chord.ChordDataLoader;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.linx.LinxDataLoader;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.protect.purple.PurpleDataLoader;
import com.hartwig.hmftools.protect.purple.ReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariantFactory;
import com.hartwig.hmftools.protect.purple.ReportableVariantSource;
import com.hartwig.hmftools.protect.virusinterpreter.VirusInterpreterData;
import com.hartwig.hmftools.protect.virusinterpreter.VirusInterpreterDataLoader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GenomicAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(GenomicAnalyzer.class);

    @NotNull
    private final GermlineReportingModel germlineReportingModel;

    public GenomicAnalyzer(@NotNull final GermlineReportingModel germlineReportingModel) {
        this.germlineReportingModel = germlineReportingModel;
    }

    @NotNull
    public GenomicAnalysis run(@NotNull String tumorSampleId, @Nullable String referenceSampleId, @NotNull PatientReporterConfig config,
            @NotNull LimsGermlineReportingLevel germlineReportingLevel) throws IOException {
        PurpleData purpleData = PurpleDataLoader.load(tumorSampleId, referenceSampleId,
                config.purpleQcFile(),
                config.purplePurityTsv(),
                config.purpleSomaticDriverCatalogTsv(),
                config.purpleSomaticVariantVcf(),
                config.purpleGermlineDriverCatalogTsv(),
                config.purpleGermlineVariantVcf(),
                null,
                config.purpleSomaticCopyNumberTsv());

        LinxData linxData = LinxDataLoader.load(config.linxFusionTsv(), config.linxBreakendTsv(), config.linxDriverCatalogTsv());
        VirusInterpreterData virusInterpreterData = VirusInterpreterDataLoader.load(config.annotatedVirusTsv());
        List<PeachGenotype> peachGenotypes = loadPeachData(config.peachGenotypeTsv());

        List<ReportableVariant> reportableVariants =
                ReportableVariantFactory.mergeVariantLists(purpleData.germlineVariants(), purpleData.somaticVariants());

        Map<ReportableVariant, Boolean> notifyGermlineStatusPerVariant =
                determineNotify(reportableVariants, germlineReportingModel, germlineReportingLevel);

        ChordAnalysis chordAnalysis = ChordDataLoader.load(config.chordPredictionTxt());

        List<ProtectEvidence> reportableEvidenceItems = extractReportableEvidenceItems(config.protectEvidenceTsv());
        List<ProtectEvidence> nonTrialsOnLabel = ReportableEvidenceItemFactory.extractNonTrialsOnLabel(reportableEvidenceItems);
        List<ProtectEvidence> trialsOnLabel = ClinicalTrialFactory.extractOnLabelTrials(reportableEvidenceItems);
        List<ProtectEvidence> nonTrialsOffLabel = ReportableEvidenceItemFactory.extractNonTrialsOffLabel(reportableEvidenceItems);

        return ImmutableGenomicAnalysis.builder()
                .impliedPurity(purpleData.purity())
                .hasReliablePurity(purpleData.hasReliablePurity())
                .hasReliableQuality(purpleData.hasReliableQuality())
                .averageTumorPloidy(purpleData.ploidy())
                .tumorSpecificEvidence(nonTrialsOnLabel)
                .clinicalTrials(trialsOnLabel)
                .offLabelEvidence(nonTrialsOffLabel)
                .reportableVariants(reportableVariants)
                .notifyGermlineStatusPerVariant(notifyGermlineStatusPerVariant)
                .microsatelliteIndelsPerMb(purpleData.microsatelliteIndelsPerMb())
                .microsatelliteStatus(purpleData.microsatelliteStatus())
                .tumorMutationalLoad(purpleData.tumorMutationalLoad())
                .tumorMutationalLoadStatus(purpleData.tumorMutationalLoadStatus())
                .tumorMutationalBurden(purpleData.tumorMutationalBurdenPerMb())
                .chordHrdValue(chordAnalysis.hrdValue())
                .chordHrdStatus(chordAnalysis.hrStatus())
                .gainsAndLosses(purpleData.reportableGainsLosses())
                .cnPerChromosome(purpleData.cnPerChromosome())
                .geneFusions(linxData.fusions())
                .geneDisruptions(linxData.geneDisruptions())
                .homozygousDisruptions(linxData.homozygousDisruptions())
                .reportableViruses(virusInterpreterData.reportableViruses())
                .peachGenotypes(peachGenotypes)
                .build();
    }

    @NotNull
    private static Map<ReportableVariant, Boolean> determineNotify(@NotNull List<ReportableVariant> reportableVariants,
            @NotNull GermlineReportingModel germlineReportingModel, @NotNull LimsGermlineReportingLevel germlineReportingLevel) {
        Map<ReportableVariant, Boolean> notifyGermlineStatusPerVariant = Maps.newHashMap();

        Set<String> germlineGenesWithIndependentHits = Sets.newHashSet();
        for (ReportableVariant variant : reportableVariants) {
            if (variant.source() == ReportableVariantSource.GERMLINE && hasOtherGermlineVariantWithDifferentPhaseSet(reportableVariants,
                    variant)) {
                germlineGenesWithIndependentHits.add(variant.gene());
            }
        }

        for (ReportableVariant variant : reportableVariants) {
            boolean notify = false;
            if (variant.source() == ReportableVariantSource.GERMLINE) {
                notify = germlineReportingModel.notifyGermlineVariant(variant, germlineReportingLevel, germlineGenesWithIndependentHits);
            }
            notifyGermlineStatusPerVariant.put(variant, notify);
        }

        return notifyGermlineStatusPerVariant;
    }

    public static boolean hasOtherGermlineVariantWithDifferentPhaseSet(@NotNull List<ReportableVariant> variants,
            @NotNull ReportableVariant variantToCompareWith) {
        Integer phaseSetToCompareWith = variantToCompareWith.localPhaseSet();
        for (ReportableVariant variant : variants) {
            if (!variant.equals(variantToCompareWith) && variant.gene().equals(variantToCompareWith.gene())
                    && variant.source() == ReportableVariantSource.GERMLINE && (phaseSetToCompareWith == null
                    || !phaseSetToCompareWith.equals(variant.localPhaseSet()))) {
                return true;
            }
        }

        return false;
    }

    @NotNull
    private static List<PeachGenotype> loadPeachData(@NotNull String peachGenotypeTsv) throws IOException {
        LOGGER.info("Loading peach genotypes from {}", new File(peachGenotypeTsv).getParent());
        List<PeachGenotype> peachGenotypes = PeachGenotypeFile.read(peachGenotypeTsv);
        LOGGER.info(" Loaded {} reportable genotypes from {}", peachGenotypes.size(), peachGenotypeTsv);
        return peachGenotypes;
    }

    @NotNull
    private static List<ProtectEvidence> extractReportableEvidenceItems(@NotNull String protectEvidenceTsv) throws IOException {
        LOGGER.info("Loading PROTECT data from {}", new File(protectEvidenceTsv).getParent());
        List<ProtectEvidence> evidences = ProtectEvidenceFile.read(protectEvidenceTsv);

        List<ProtectEvidence> reportableEvidenceItems = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (evidence.reported()) {
                reportableEvidenceItems.add(evidence);
            }
        }
        LOGGER.info(" Loaded {} reportable evidence items from {}", reportableEvidenceItems.size(), protectEvidenceTsv);
        return reportableEvidenceItems;
    }
}
