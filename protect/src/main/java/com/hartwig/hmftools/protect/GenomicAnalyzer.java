package com.hartwig.hmftools.protect;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordFileReader;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.clinical.PatientTumorLocation;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
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
import com.hartwig.hmftools.protect.structural.SvAnalysis;
import com.hartwig.hmftools.protect.structural.SvAnalyzer;
import com.hartwig.hmftools.protect.variants.ReportableVariantAnalysis;
import com.hartwig.hmftools.protect.variants.ReportableVariantAnalyzer;
import com.hartwig.hmftools.protect.variants.germline.DriverGermlineVariant;
import com.hartwig.hmftools.protect.variants.germline.FilterGermlineVariants;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.protect.variants.somatic.DriverSomaticVariant;
import com.hartwig.hmftools.protect.variants.somatic.SomaticVariantAnalyzer;
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
    public GenomicAnalysis run(@NotNull String tumorSampleId, @Nullable PatientTumorLocation patientTumorLocation,
            @NotNull LimsGermlineReportingLevel germlineReportingLevel, boolean reportViralInsertions, @NotNull String purplePurityTsv,
            @NotNull String purpleQCFile, @NotNull String purpleGeneCnvTsv, @NotNull String purpleDriverCatalogTsv,
            @NotNull String purpleSomaticVariantVcf, @NotNull String bachelorTsv, @NotNull String linxFusionTsv,
            @NotNull String linxBreakendTsv, @NotNull String linxViralInsertionTsv, @NotNull String linxDriversTsv,
            @NotNull String chordPredictionTxt) throws IOException {
        List<DriverCatalog> purpleDriverCatalog = readDriverCatalog(purpleDriverCatalogTsv);
        PurpleAnalysis purpleAnalysis =
                analyzePurple(purplePurityTsv, purpleQCFile, purpleGeneCnvTsv, patientTumorLocation, purpleDriverCatalog);
        List<DriverSomaticVariant> driverSomaticVariants =
                analyzeSomaticVariants(tumorSampleId, purpleSomaticVariantVcf, purpleDriverCatalog);

        ChordAnalysis chordAnalysis = analyzeChord(chordPredictionTxt);
        ChordStatus chordStatus = chordAnalysis.hrStatus();

        List<DriverGermlineVariant> driverGermlineVariants =
                analyzeGermlineVariants(bachelorTsv, purpleAnalysis, driverSomaticVariants, chordStatus, germlineReportingLevel);

        ReportableVariantAnalysis reportableVariantsAnalysis = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(
                driverSomaticVariants,
                driverGermlineVariants,
                germlineReportingModel,
                germlineReportingLevel,
                actionabilityAnalyzer,
                patientTumorLocation);

        SvAnalysis svAnalysis = analyzeStructuralVariants(linxFusionTsv, linxBreakendTsv, patientTumorLocation);
        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = extractHomozygousDisruptionsFromLinxDrivers(linxDriversTsv);
        List<ViralInsertion> viralInsertions = analyzeViralInsertions(linxViralInsertionTsv, reportViralInsertions);

        List<EvidenceItem> allEvidenceItems = Lists.newArrayList();
        allEvidenceItems.addAll(reportableVariantsAnalysis.evidenceItems());
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
                .reportableVariants(reportableVariantsAnalysis.variantsToReport())
                .microsatelliteIndelsPerMb(purpleAnalysis.purpleSignatures().microsatelliteIndelsPerMb())
                .microsatelliteStatus(purpleAnalysis.purpleSignatures().microsatelliteStatus())
                .tumorMutationalLoad(purpleAnalysis.purpleSignatures().tumorMutationalLoad())
                .tumorMutationalLoadStatus(purpleAnalysis.purpleSignatures().tumorMutationalLoadStatus())
                .tumorMutationalBurden(purpleAnalysis.purpleSignatures().tumorMutationalBurdenPerMb())
                .chordHrdValue(chordAnalysis.hrdValue())
                .chordHrdStatus(chordStatus)
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
    private PurpleAnalysis analyzePurple(@NotNull String purplePurityTsv, @NotNull String purpleQCFile, @NotNull String purpleGeneCnvTsv,
            @Nullable PatientTumorLocation patientTumorLocation, @NotNull List<DriverCatalog> purpleDriverCatalog) throws IOException {
        PurityContext purityContext = PurityContextFile.readWithQC(purpleQCFile, purplePurityTsv);
        LOGGER.info("Loaded purple sample data from {}", purplePurityTsv);
        LOGGER.info(" Purple purity: {}", new DecimalFormat("#'%'").format(purityContext.bestFit().purity() * 100));
        LOGGER.info(" Purple average tumor ploidy: {}", purityContext.bestFit().ploidy());
        LOGGER.info(" Purple fit method: {}", purityContext.method());
        LOGGER.info(" WGD happened: {}", purityContext.wholeGenomeDuplication() ? "yes" : "no");

        PurpleQC purpleQC = purityContext.qc();
        LOGGER.info("Loaded purple QC data from {}", purpleQCFile);
        LOGGER.info(" Purple QC status: {}", purpleQC.toString());

        List<GeneCopyNumber> exomeGeneCopyNumbers = GeneCopyNumberFile.read(purpleGeneCnvTsv);
        LOGGER.info("Loaded {} gene copy numbers from {}", exomeGeneCopyNumbers.size(), purpleGeneCnvTsv);

        return PurpleAnalyzer.run(purityContext,
                purpleQC,
                exomeGeneCopyNumbers,
                actionabilityAnalyzer,
                patientTumorLocation,
                purpleDriverCatalog);
    }

    @NotNull
    private static List<DriverSomaticVariant> analyzeSomaticVariants(@NotNull String sample, @NotNull String somaticVariantVcf,
            @NotNull List<DriverCatalog> purpleDriverCatalog) throws IOException {
        List<SomaticVariant> variants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(sample, somaticVariantVcf);
        LOGGER.info("Loaded {} PASS somatic variants from {}", variants.size(), somaticVariantVcf);

        return SomaticVariantAnalyzer.run(variants, purpleDriverCatalog);
    }

    @NotNull
    private List<DriverGermlineVariant> analyzeGermlineVariants(@NotNull String bachelorTsv, @NotNull PurpleAnalysis purpleAnalysis,
            @NotNull List<DriverSomaticVariant> driverSomaticVariants, @NotNull ChordStatus chordStatus,
            @NotNull LimsGermlineReportingLevel germlineChoice) throws IOException {
        List<ReportableGermlineVariant> variants = ReportableGermlineVariantFile.read(bachelorTsv);

        LOGGER.info("Loaded {} reportable germline variants from {}", variants.size(), bachelorTsv);

        if (germlineChoice != LimsGermlineReportingLevel.NO_REPORTING) {
            LOGGER.info(" Patient has given the following germline consent: '{}'", germlineChoice);
            return FilterGermlineVariants.filterGermlineVariantsForReporting(variants,
                    germlineReportingModel,
                    purpleAnalysis.exomeGeneCopyNumbers(),
                    driverSomaticVariants,
                    chordStatus);
        } else {
            LOGGER.info(" No consent has been given for germline reporting. No germline variants will be reported!");
            return Lists.newArrayList();
        }
    }

    @NotNull
    private SvAnalysis analyzeStructuralVariants(@NotNull String linxFusionTsv, @NotNull String linxBreakendTsv,
            @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        List<LinxFusion> linxFusions = LinxFusion.read(linxFusionTsv);
        List<LinxFusion> linxReportedFusions = linxFusions.stream().filter(x -> x.reported()).collect(Collectors.toList());

        LOGGER.info("Loaded {} fusions from {}", linxReportedFusions.size(), linxFusionTsv);

        List<LinxBreakend> linxBreakends = LinxBreakend.read(linxBreakendTsv);
        List<LinxBreakend> linxReportedBreakends = linxBreakends.stream().filter(x -> x.reportedDisruption()).collect(Collectors.toList());

        LOGGER.info("Loaded {} disruptions from {}", linxReportedBreakends.size(), linxBreakendTsv);

        return SvAnalyzer.run(linxReportedFusions, linxReportedBreakends, actionabilityAnalyzer, patientTumorLocation);
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
    private static List<ViralInsertion> analyzeViralInsertions(@NotNull String linxViralInsertionTsv, boolean reportViralInsertions)
            throws IOException {
        List<LinxViralInsertion> viralInsertionList = LinxViralInsertion.read(linxViralInsertionTsv);
        LOGGER.info("Loaded {} viral insertions from {}", viralInsertionList.size(), linxViralInsertionTsv);

        if (reportViralInsertions) {
            List<ViralInsertion> reportableViralInsertions = ViralInsertionAnalyzer.analyzeViralInsertions(viralInsertionList);
            LOGGER.info(" Patient has given consent for viral insertion reporting. Found {} reportable viral insertions.",
                    reportableViralInsertions.size());
            return reportableViralInsertions;
        } else {
            LOGGER.info(" No consent has been given for viral insertions. No viral insertions will be reported!");
            return null;
        }
    }
}
