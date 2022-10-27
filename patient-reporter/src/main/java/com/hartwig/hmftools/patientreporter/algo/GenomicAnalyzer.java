package com.hartwig.hmftools.patientreporter.algo;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordDataFile;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.hla.HlaAllelesReportingData;
import com.hartwig.hmftools.common.hla.HlaAllelesReportingFactory;
import com.hartwig.hmftools.common.hla.LilacSummaryData;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.linx.GeneDisruption;
import com.hartwig.hmftools.common.linx.GeneDisruptionFactory;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.linx.LinxDataLoader;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceFile;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.loader.PurpleData;
import com.hartwig.hmftools.common.purple.loader.PurpleDataLoader;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantFactory;
import com.hartwig.hmftools.common.variant.ReportableVariantSource;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.common.virus.VirusInterpreterDataLoader;
import com.hartwig.hmftools.patientreporter.PatientReporterConfig;
import com.hartwig.hmftools.patientreporter.actionability.ClinicalTrialFactory;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.algo.orange.BreakendSelector;
import com.hartwig.hmftools.patientreporter.algo.orange.LossOfHeterozygositySelector;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingModel;

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
            @NotNull LimsGermlineReportingLevel germlineReportingLevel, @NotNull KnownFusionCache knownFusionCache) throws IOException {
        PurpleData purpleData = PurpleDataLoader.load(tumorSampleId,
                referenceSampleId,
                null,
                config.purpleQcFile(),
                config.purplePurityTsv(),
                config.purpleSomaticDriverCatalogTsv(),
                config.purpleSomaticVariantVcf(),
                config.purpleGermlineDriverCatalogTsv(),
                config.purpleGermlineVariantVcf(),
                null,
                config.purpleSomaticCopyNumberTsv(),
                null,
                config.refGenomeVersion());

        LinxData linxData = LinxDataLoader.load(config.linxSvsTsv(),
                config.linxFusionTsv(),
                config.linxBreakendTsv(),
                config.linxDriverCatalogTsv(),
                null,
                null);

        List<GeneDisruption> additionalSuspectBreakends =
                GeneDisruptionFactory.convert(BreakendSelector.selectInterestingUnreportedBreakends(linxData.allBreakends(),
                        linxData.reportableFusions(),
                        knownFusionCache), linxData.allStructuralVariants());

        List<GeneDisruption> reportableGeneDisruptions = Lists.newArrayList();
        reportableGeneDisruptions.addAll(linxData.reportableGeneDisruptions());
        reportableGeneDisruptions.addAll(additionalSuspectBreakends);

        ChordData chordAnalysis = ChordDataFile.read(config.chordPredictionTxt(), true);

        List<GeneCopyNumber> suspectGeneCopyNumbersHRDWithLOH =
                LossOfHeterozygositySelector.selectHRDGenesWithLOH(purpleData.allSomaticGeneCopyNumbers(), chordAnalysis.hrStatus());
        LOGGER.info(" Found an additional {} suspect gene copy numbers HRD with LOH", suspectGeneCopyNumbersHRDWithLOH.size());

        List<GeneCopyNumber> suspectGeneCopyNumbersMSIWithLOH =
                LossOfHeterozygositySelector.selectMSIGenesWithLOH(purpleData.allSomaticGeneCopyNumbers(),
                        purpleData.microsatelliteStatus());
        LOGGER.info(" Found an additional {} suspect gene copy numbers HRD with LOH", suspectGeneCopyNumbersMSIWithLOH.size());

        VirusInterpreterData virusInterpreterData = VirusInterpreterDataLoader.load(config.annotatedVirusTsv());

        List<ReportableVariant> reportableVariants =
                ReportableVariantFactory.mergeVariantLists(purpleData.reportableGermlineVariants(), purpleData.reportableSomaticVariants());

        Map<ReportableVariant, Boolean> notifyGermlineStatusPerVariant =
                determineNotify(reportableVariants, germlineReportingModel, germlineReportingLevel);

        LilacSummaryData lilacSummaryData = LilacSummaryData.load(config.lilacQcCsv(), config.lilacResultCsv());
        HlaAllelesReportingData hlaReportingData = HlaAllelesReportingFactory.convertToReportData(lilacSummaryData, purpleData.hasReliablePurity());

        List<ProtectEvidence> reportableEvidenceItems = extractReportableEvidenceItems(config.protectEvidenceTsv());
        List<ProtectEvidence> nonTrialsOnLabel = ReportableEvidenceItemFactory.extractNonTrialsOnLabel(reportableEvidenceItems);
        List<ProtectEvidence> trialsOnLabel = ClinicalTrialFactory.extractOnLabelTrials(reportableEvidenceItems);
        List<ProtectEvidence> nonTrialsOffLabel = ReportableEvidenceItemFactory.extractNonTrialsOffLabel(reportableEvidenceItems);

        return ImmutableGenomicAnalysis.builder()
                .purpleQCStatus(purpleData.qc().status())
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
                .hrdValue(chordAnalysis.hrdValue())
                .hrdStatus(chordAnalysis.hrStatus())
                .gainsAndLosses(purpleData.reportableSomaticGainsLosses())
                .cnPerChromosome(purpleData.copyNumberPerChromosome())
                .geneFusions(linxData.reportableFusions())
                .geneDisruptions(reportableGeneDisruptions)
                .homozygousDisruptions(linxData.homozygousDisruptions())
                .reportableViruses(virusInterpreterData.reportableViruses())
                .hlaAlleles(hlaReportingData)
                .suspectGeneCopyNumbersHRDWithLOH(suspectGeneCopyNumbersHRDWithLOH)
                .suspectGeneCopyNumbersMSIWithLOH(suspectGeneCopyNumbersMSIWithLOH)
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
