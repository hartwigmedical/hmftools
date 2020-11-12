package com.hartwig.hmftools.protect;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordFileReader;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariantFile;
import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertion;
import com.hartwig.hmftools.protect.actionability.ClinicalTrialFactory;
import com.hartwig.hmftools.protect.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.protect.homozygousdisruption.HomozygousDisruptionAnalyzer;
import com.hartwig.hmftools.protect.homozygousdisruption.ReportableHomozygousDisruption;
import com.hartwig.hmftools.protect.purple.PurpleAnalysis;
import com.hartwig.hmftools.protect.purple.PurpleAnalyzer;
import com.hartwig.hmftools.protect.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.protect.structural.SvAnalysis;
import com.hartwig.hmftools.protect.structural.SvAnalyzer;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantAnalysis;
import com.hartwig.hmftools.protect.variants.ReportableVariantAnalyzer;
import com.hartwig.hmftools.protect.variants.ReportableVariantFactory;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.protect.viralinsertion.ViralInsertion;
import com.hartwig.hmftools.protect.viralinsertion.ViralInsertionAnalyzer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GenomicAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(GenomicAnalyzer.class);

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
        List<DriverCatalog> purpleDriverCatalog = readDriverCatalog(purpleDriverCatalogTsv);
        PurpleAnalysis purpleAnalysis = analyzePurple(purplePurityTsv, purpleQCFile, patientPrimaryTumor, purpleDriverCatalog);
        List<ReportableVariant> reportableSomaticVariants =
                analyzeSomaticVariants(tumorSampleId, purpleSomaticVariantVcf, purpleDriverCatalog);

        SvAnalysis svAnalysis = analyzeStructuralVariants(linxFusionTsv, linxBreakendTsv, patientPrimaryTumor);
        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = extractHomozygousDisruptionsFromLinxDrivers(linxDriversTsv);
        List<ViralInsertion> viralInsertions = analyzeViralInsertions(linxViralInsertionTsv);

        List<ReportableVariant> reportableGermlineVariants = analyzeGermlineVariants(bachelorTsv,
                reportableSomaticVariants,
                purpleAnalysis.reportableGainsAndLosses(),
                reportableHomozygousDisruptions,
                svAnalysis.reportableDisruptions());

        ReportableVariantAnalysis reportableVariantAnalysis = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(
                reportableSomaticVariants,
                reportableGermlineVariants,
                actionabilityAnalyzer, patientPrimaryTumor);

        ChordAnalysis chordAnalysis = analyzeChord(chordPredictionTxt);

        List<EvidenceItem> allEvidenceItems = Lists.newArrayList();
        allEvidenceItems.addAll(reportableVariantAnalysis.evidenceItems());
        allEvidenceItems.addAll(purpleAnalysis.evidenceItems());
        allEvidenceItems.addAll(svAnalysis.evidenceItems());

        List<EvidenceItem> allEvidenceItemsFiltered = ReportableEvidenceItemFactory.filterBlacklistedEvidence(allEvidenceItems);
        List<EvidenceItem> nonTrials = ReportableEvidenceItemFactory.extractNonTrials(allEvidenceItemsFiltered);

        return ImmutableGenomicAnalysis.builder()
                .impliedPurity(purpleAnalysis.purity())
                .hasReliablePurity(purpleAnalysis.hasReliablePurity())
                .hasReliableQuality(purpleAnalysis.hasReliableQuality())
                .averageTumorPloidy(purpleAnalysis.ploidy())
                .tumorSpecificEvidence(nonTrials.stream().filter(EvidenceItem::isOnLabel).collect(Collectors.toList()))
                .clinicalTrials(ClinicalTrialFactory.extractOnLabelTrials(allEvidenceItemsFiltered))
                .offLabelEvidence(nonTrials.stream().filter(item -> !item.isOnLabel()).collect(Collectors.toList()))
                .reportableVariants(reportableVariantAnalysis.variantsToReport())
                .microsatelliteIndelsPerMb(purpleAnalysis.purpleSignatures().microsatelliteIndelsPerMb())
                .microsatelliteStatus(purpleAnalysis.purpleSignatures().microsatelliteStatus())
                .tumorMutationalLoad(purpleAnalysis.purpleSignatures().tumorMutationalLoad())
                .tumorMutationalLoadStatus(purpleAnalysis.purpleSignatures().tumorMutationalLoadStatus())
                .tumorMutationalBurden(purpleAnalysis.purpleSignatures().tumorMutationalBurdenPerMb())
                .chordHrdValue(chordAnalysis.hrdValue())
                .chordHrdStatus(chordAnalysis.hrStatus())
                .gainsAndLosses(purpleAnalysis.reportableGainsAndLosses())
                .geneFusions(svAnalysis.reportableFusions())
                .geneDisruptions(svAnalysis.reportableDisruptions())
                .homozygousDisruptions(reportableHomozygousDisruptions)
                .viralInsertions(viralInsertions)
                .build();
    }

    @NotNull
    private static List<DriverCatalog> readDriverCatalog(@NotNull String purpleDriverCatalogTsv) throws IOException {
        List<DriverCatalog> driverCatalog = DriverCatalogFile.read(purpleDriverCatalogTsv);
        LOGGER.info("Loaded {} purple driver catalog records from {}", driverCatalog.size(), purpleDriverCatalogTsv);
        return driverCatalog;
    }

    @NotNull
    private PurpleAnalysis analyzePurple(@NotNull String purplePurityTsv, @NotNull String purpleQCFile,
            @Nullable PatientPrimaryTumor patientPrimaryTumor, @NotNull List<DriverCatalog> purpleDriverCatalog) throws IOException {
        PurityContext purityContext = PurityContextFile.readWithQC(purpleQCFile, purplePurityTsv);
        LOGGER.info("Loaded purple sample data from {}", purplePurityTsv);
        LOGGER.info(" Purple purity: {}", new DecimalFormat("#'%'").format(purityContext.bestFit().purity() * 100));
        LOGGER.info(" Purple average tumor ploidy: {}", purityContext.bestFit().ploidy());
        LOGGER.info(" Purple fit method: {}", purityContext.method());
        LOGGER.info(" WGD happened: {}", purityContext.wholeGenomeDuplication() ? "yes" : "no");

        PurpleQC purpleQC = purityContext.qc();
        LOGGER.info("Loaded purple QC data from {}", purpleQCFile);
        LOGGER.info(" Purple QC status: {}", purpleQC.toString());

        return PurpleAnalyzer.run(purityContext, purpleQC, actionabilityAnalyzer, patientPrimaryTumor, purpleDriverCatalog);
    }

    @NotNull
    private static List<ReportableVariant> analyzeSomaticVariants(@NotNull String sample, @NotNull String somaticVariantVcf,
            @NotNull List<DriverCatalog> purpleDriverCatalog) throws IOException {
        List<SomaticVariant> variants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(sample, somaticVariantVcf);
        LOGGER.info("Loaded {} PASS somatic variants from {}", variants.size(), somaticVariantVcf);

        return ReportableVariantFactory.reportableSomaticVariants(variants, purpleDriverCatalog);
    }

    @NotNull
    private List<ReportableVariant> analyzeGermlineVariants(@NotNull String bachelorTsv,
            @NotNull List<ReportableVariant> reportableSomaticVariants, @NotNull List<ReportableGainLoss> reportableGainLosses,
            @NotNull List<ReportableHomozygousDisruption> reportableHomozygousDisruptions,
            @NotNull List<ReportableGeneDisruption> reportableGeneDisruptions) throws IOException {
        List<ReportableGermlineVariant> germlineVariants = ReportableGermlineVariantFile.read(bachelorTsv);

        LOGGER.info("Loaded {} reportable germline variants from {}", germlineVariants.size(), bachelorTsv);

        Set<String> somaticGenes = Sets.newHashSet();
        reportableHomozygousDisruptions.stream().map(ReportableHomozygousDisruption::gene).forEach(somaticGenes::add);
        reportableGeneDisruptions.stream().map(ReportableGeneDisruption::gene).forEach(somaticGenes::add);
        reportableSomaticVariants.stream().map(ReportableVariant::gene).forEach(somaticGenes::add);
        reportableGainLosses.stream()
                .filter(x -> !x.interpretation().equals(CopyNumberInterpretation.GAIN))
                .map(ReportableGainLoss::gene)
                .forEach(somaticGenes::add);

        List<ReportableVariant> reportableVariants =
                ReportableVariantFactory.reportableGermlineVariants(germlineVariants, somaticGenes, germlineReportingModel);

        LOGGER.info(" Filtered to {} driver germline variants", reportableVariants.size());
        return reportableVariants;
    }

    @NotNull
    private SvAnalysis analyzeStructuralVariants(@NotNull String linxFusionTsv, @NotNull String linxBreakendTsv,
            @Nullable PatientPrimaryTumor patientPrimaryTumor) throws IOException {
        List<LinxFusion> linxFusions = LinxFusion.read(linxFusionTsv);
        List<LinxFusion> linxReportedFusions = linxFusions.stream().filter(x -> x.reported()).collect(Collectors.toList());

        LOGGER.info("Loaded {} fusions from {}", linxReportedFusions.size(), linxFusionTsv);

        List<LinxBreakend> linxBreakends = LinxBreakend.read(linxBreakendTsv);
        List<LinxBreakend> linxReportedBreakends = linxBreakends.stream().filter(x -> x.reportedDisruption()).collect(Collectors.toList());

        LOGGER.info("Loaded {} disruptions from {}", linxReportedBreakends.size(), linxBreakendTsv);

        return SvAnalyzer.run(linxReportedFusions, linxReportedBreakends, actionabilityAnalyzer, patientPrimaryTumor);
    }

    @NotNull
    private static ChordAnalysis analyzeChord(@NotNull String chordPredictionTxt) throws IOException {
        ChordAnalysis chord = ChordFileReader.read(chordPredictionTxt);
        LOGGER.info("Loaded CHORD analysis from {}", chordPredictionTxt);
        return chord;
    }

    @NotNull
    private static List<ReportableHomozygousDisruption> extractHomozygousDisruptionsFromLinxDrivers(@NotNull String linxDriversTsv)
            throws IOException {
        return HomozygousDisruptionAnalyzer.extractFromLinxDriversTsv(linxDriversTsv);
    }

    @Nullable
    private static List<ViralInsertion> analyzeViralInsertions(@NotNull String linxViralInsertionTsv) throws IOException {
        List<LinxViralInsertion> viralInsertionList = LinxViralInsertion.read(linxViralInsertionTsv);
        LOGGER.info("Loaded {} viral insertions from {}", viralInsertionList.size(), linxViralInsertionTsv);

        List<ViralInsertion> reportableViralInsertions = ViralInsertionAnalyzer.analyzeViralInsertions(viralInsertionList);
        LOGGER.info(" Viral insertions have been consolidated into {} reportable viral insertions.", reportableViralInsertions.size());

        return reportableViralInsertions;
    }
}
