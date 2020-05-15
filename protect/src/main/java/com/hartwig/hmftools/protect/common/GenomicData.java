package com.hartwig.hmftools.protect.common;

import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusionFile;
import com.hartwig.hmftools.protect.report.chord.ChordAnalysis;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GenomicData {

    private static final Logger LOGGER = LogManager.getLogger(GenomicData.class);

    private GenomicData() {

    }

    public static double extractPloidy(@NotNull String purplePurityTsv) throws IOException {
        LOGGER.info("Reading purple purity from {}", purplePurityTsv);
        PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);
        double ploidy = purityContext.bestFit().ploidy();
        LOGGER.info(" Sample ploidy: {}", ploidy);
        return ploidy;
    }

    @NotNull
    public static List<? extends Variant> readPassSomaticVariants(@NotNull String sampleId, @NotNull String somaticVariantVcf)
            throws IOException {
        LOGGER.info("Reading somatic variants from {}", somaticVariantVcf);
        List<? extends Variant> passSomaticVariants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(sampleId, somaticVariantVcf);
        LOGGER.info(" Loaded {} PASS somatic variants", passSomaticVariants.size());
        return passSomaticVariants;
    }

    @NotNull
    public static List<GeneCopyNumber> readGeneCopyNumbers(@NotNull String purpleGeneCnvTsv) throws IOException {
        LOGGER.info("Reading gene copy numbers from {}", purpleGeneCnvTsv);
        List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(purpleGeneCnvTsv);
        LOGGER.info(" Loaded {} gene copy numbers", geneCopyNumbers.size());
        return geneCopyNumbers;
    }

    @NotNull
    public static List<ReportableGeneFusion> readGeneFusions(@NotNull String linxFusionTsv) throws IOException {
        LOGGER.info("Reading gene fusions from {}", linxFusionTsv);
        List<ReportableGeneFusion> fusions = ReportableGeneFusionFile.read(linxFusionTsv);
        LOGGER.info(" Loaded {} fusions", fusions.size());
        return fusions;
    }

    @NotNull
    public static List<ReportableGermlineVariant> analyzeGermlineVariants(@NotNull String sampleBarcode, @NotNull String bachelorTsv,
            @NotNull CopyNumberAnalysis copyNumberAnalysis, @NotNull SomaticVariantAnalysis somaticVariantAnalysis,
            @NotNull ChordAnalysis chordAnalysis, @NotNull DriverGeneView driverGeneView, @NotNull Lims lims,
            @NotNull GermlineReportingModel germlineReportingModel) throws IOException {

        List<GermlineVariant> variants =
                BachelorFile.loadBachelorTsv(bachelorTsv).stream().filter(GermlineVariant::passFilter).collect(Collectors.toList());
        LOGGER.info("Loaded {} PASS germline variants from {}", variants.size(), bachelorTsv);

        LimsGermlineReportingLevel germlineChoice = lims.germlineReportingChoice(sampleBarcode);
        if (germlineChoice == LimsGermlineReportingLevel.NO_REPORTING) {
            LOGGER.info(" No germline reporting choice known. No germline variants will be reported!");
            return Lists.newArrayList();
        } else {
            LOGGER.info(" Patient has given the following germline consent: {}", germlineChoice);
            return filterGermlineVariantsForReporting(variants,
                    driverGeneView,
                    germlineReportingModel,
                    copyNumberAnalysis.exomeGeneCopyNumbers(),
                    somaticVariantAnalysis.variantsToReport(),
                    chordAnalysis);
        }
    }

    @NotNull
    public static List<ReportableGermlineVariant> filterGermlineVariantsForReporting(@NotNull List<GermlineVariant> germlineVariants,
            @NotNull DriverGeneView driverGeneView, @NotNull GermlineReportingModel germlineReportingModel,
            @NotNull List<GeneCopyNumber> allGeneCopyNumbers, @NotNull List<SomaticVariant> variantsToReport,
            @NotNull ChordAnalysis chordAnalysis) {
        List<ReportableGermlineVariant> reportableGermlineVariants = Lists.newArrayList();

        Set<String> reportableGermlineGenes = germlineReportingModel.reportableGermlineGenes();
        for (GermlineVariant germlineVariant : germlineVariants) {
            assert germlineVariant.passFilter();

            if (reportableGermlineGenes.contains(germlineVariant.gene())) {
                // Note: Reportable germline genes may not necessarily be present in driverGeneView!
                if (driverGeneView.category(germlineVariant.gene()) == DriverCategory.ONCO) {
                    // Report all germline variants on reportable oncogenes.
                    reportableGermlineVariants.add(reportableGermlineVariantWithDriverLikelihood(germlineVariant, 1.0));
                } else {
                    // Only report germline variants on TSGs if there is a 2nd hit or CHORD suggests HRD
                    boolean filterBiallelic = germlineVariant.biallelic();

                    boolean filterMinCopyNumberTumor = false;
                    GeneCopyNumber geneCopyNumber = lookupGeneCopyNumber(allGeneCopyNumbers, germlineVariant.gene());
                    if (Math.round(geneCopyNumber.minCopyNumber()) <= 1 && (Math.round(germlineVariant.adjustedCopyNumber()) >= 2)) {
                        filterMinCopyNumberTumor = true;
                    }

                    boolean filterSomaticVariantInSameGene = false;
                    for (SomaticVariant variant : variantsToReport) {
                        if (variant.gene().equals(germlineVariant.gene())) {
                            filterSomaticVariantInSameGene = true;
                        }
                    }

                    boolean filterGermlineVariantInSameGene = false;
                    for (GermlineVariant variant : germlineVariants) {
                        if (variant != germlineVariant && variant.gene().equals(germlineVariant.gene())) {
                            filterGermlineVariantInSameGene = true;
                        }
                    }

                    if (filterBiallelic || filterSomaticVariantInSameGene || filterGermlineVariantInSameGene) {
                        reportableGermlineVariants.add(reportableGermlineVariantWithDriverLikelihood(germlineVariant, 1.0));
                    } else if (filterMinCopyNumberTumor || chordAnalysis.predictedResponseValue()) {
                        reportableGermlineVariants.add(reportableGermlineVariantWithDriverLikelihood(germlineVariant, 0.5));
                    }
                }
            }
        }
        return reportableGermlineVariants;
    }

    @NotNull
    private static ReportableGermlineVariant reportableGermlineVariantWithDriverLikelihood(@NotNull GermlineVariant germlineVariant,
            double driverLikelihood) {
        return ImmutableReportableGermlineVariant.builder().variant(germlineVariant).driverLikelihood(driverLikelihood).build();
    }

    @NotNull
    private static GeneCopyNumber lookupGeneCopyNumber(@NotNull List<GeneCopyNumber> allGeneCopyNumbers, @NotNull String gene) {
        for (GeneCopyNumber geneCopyNumber : allGeneCopyNumbers) {
            if (geneCopyNumber.gene().equals(gene)) {
                return geneCopyNumber;
            }
        }

        throw new IllegalStateException("Could not find gene copy number for gene: " + gene);
    }

    @NotNull
    public static ReportableVariantAnalysis somaticAndGermlineVariantsTogether(@NotNull List<SomaticVariant> somaticVariantsReport,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull DriverGeneView driverGeneView,
            @NotNull List<ReportableGermlineVariant> germlineVariantsToReport, @NotNull GermlineReportingModel germlineReportingModel,
            @NotNull LimsGermlineReportingLevel germlineReportingChoice) {
        List<ReportableVariant> allReportableVariants = mergeSomaticAndGermlineVariants(somaticVariantsReport,
                driverCatalog,
                driverGeneView,
                germlineVariantsToReport,
                germlineReportingModel,
                germlineReportingChoice);

        return ImmutableReportableVariantAnalysis.builder().variantsToReport(allReportableVariants).build();

    }

    @NotNull
    private static List<ReportableVariant> mergeSomaticAndGermlineVariants(@NotNull List<SomaticVariant> somaticVariantsReport,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull DriverGeneView driverGeneView,
            @NotNull List<ReportableGermlineVariant> germlineVariantsToReport, @NotNull GermlineReportingModel germlineReportingModel,
            @NotNull LimsGermlineReportingLevel germlineReportingChoice) {
        List<ReportableVariant> allReportableVariants = Lists.newArrayList();
        for (SomaticVariant variant : somaticVariantsReport) {
            DriverCategory category = driverGeneView.category(variant.gene());
            assert category != null;

            DriverCatalog catalog = catalogEntryForVariant(driverCatalog, variant.gene());
            Double driverLikelihood = null;
            if (catalog != null) {
                driverLikelihood = catalog.driverLikelihood();
                for (ReportableGermlineVariant germlineVariant : germlineVariantsToReport) {
                    if (germlineVariant.variant().gene().equals(variant.gene())) {
                        driverLikelihood = Math.max(driverLikelihood, germlineVariant.driverLikelihood());
                    }
                }
            }

            allReportableVariants.add(fromSomaticVariant(variant).driverCategory(category)
                    .driverLikelihood(driverLikelihood)
                    .notifyClinicalGeneticist(false)
                    .build());
        }

        boolean wantsToBeNotified = germlineReportingChoice == LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION;
        for (ReportableGermlineVariant germlineVariant : germlineVariantsToReport) {
            DriverCategory category = driverGeneView.category(germlineVariant.variant().gene());
            DriverCatalog catalog = catalogEntryForVariant(driverCatalog, germlineVariant.variant().gene());
            double driverLikelihood = germlineVariant.driverLikelihood();
            if (catalog != null) {
                driverLikelihood = Math.max(driverLikelihood, catalog.driverLikelihood());
            }
            allReportableVariants.add(fromGermlineVariant(germlineVariant.variant()).driverCategory(category)
                    .driverLikelihood(driverLikelihood)
                    .notifyClinicalGeneticist(wantsToBeNotified && germlineReportingModel.notifyAboutGene(germlineVariant.variant().gene()))
                    .build());
        }

        return allReportableVariants;
    }

    @Nullable
    private static DriverCatalog catalogEntryForVariant(@NotNull List<DriverCatalog> driverCatalogList, @NotNull String gene) {
        for (DriverCatalog entry : driverCatalogList) {
            if (entry.gene().equals(gene)) {
                return entry;
            }
        }
        return null;
    }

    @NotNull
    private static ImmutableReportableVariant.Builder fromGermlineVariant(@NotNull GermlineVariant variant) {
        return ImmutableReportableVariant.builder()
                .gene(variant.gene())
                .position(variant.position())
                .chromosome(variant.chromosome())
                .ref(variant.ref())
                .alt(variant.alt())
                .canonicalCodingEffect(variant.codingEffect())
                .canonicalHgvsCodingImpact(variant.hgvsCodingImpact())
                .canonicalHgvsProteinImpact(variant.hgvsProteinImpact())
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .gDNA(toGDNA(variant))
                .totalPloidy(variant.adjustedCopyNumber())
                .allelePloidy(calcAllelePloidy(variant.adjustedCopyNumber(), variant.adjustedVAF()))
                .hotspot(Hotspot.NON_HOTSPOT)
                .clonalLikelihood(1D)
                .biallelic(variant.biallelic());
    }

    @NotNull
    private static ImmutableReportableVariant.Builder fromSomaticVariant(@NotNull SomaticVariant variant) {
        return ImmutableReportableVariant.builder()
                .gene(variant.gene())
                .position(variant.position())
                .chromosome(variant.chromosome())
                .ref(variant.ref())
                .alt(variant.alt())
                .canonicalCodingEffect(variant.canonicalCodingEffect())
                .canonicalHgvsCodingImpact(variant.canonicalHgvsCodingImpact())
                .canonicalHgvsProteinImpact(variant.canonicalHgvsProteinImpact())
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .gDNA(toGDNA(variant))
                .totalPloidy(variant.adjustedCopyNumber())
                .allelePloidy(calcAllelePloidy(variant.adjustedCopyNumber(), variant.adjustedVAF()))
                .hotspot(variant.hotspot())
                .clonalLikelihood(variant.clonalLikelihood())
                .biallelic(variant.biallelic());
    }

    @NotNull
    private static String toGDNA(@NotNull GenomePosition genomePosition) {
        return genomePosition.chromosome() + ":" + genomePosition.position();
    }

    private static double calcAllelePloidy(double adjustedCopyNumber, double adjustedVAF) {
        return adjustedCopyNumber * Math.max(0, Math.min(1, adjustedVAF));
    }
}
