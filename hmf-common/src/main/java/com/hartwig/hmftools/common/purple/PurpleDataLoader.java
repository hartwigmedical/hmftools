package com.hartwig.hmftools.common.purple;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.drivercatalog.CNADrivers;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGenePanel;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.cnchromosome.CnPerChromosomeArmData;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.GenerateCnPerChromosome;
import com.hartwig.hmftools.common.purple.copynumber.ImmutableReportableGainLoss;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class PurpleDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(PurpleDataLoader.class);

    private PurpleDataLoader() {
    }

    @NotNull
    public static PurpleData load(@NotNull String tumorSample, @Nullable String referenceSample, @NotNull String qcFile,
            @NotNull String purityTsv, @NotNull String somaticDriverCatalogTsv, @NotNull String somaticVariantVcf,
            @NotNull String germlineDriverCatalogTsv, @NotNull String germlineVariantVcf, @Nullable String purpleGeneCopyNumberTsv)
            throws IOException {
        return load(tumorSample,
                referenceSample,
                qcFile,
                purityTsv,
                somaticDriverCatalogTsv,
                somaticVariantVcf,
                germlineDriverCatalogTsv,
                germlineVariantVcf,
                purpleGeneCopyNumberTsv,
                null,
                null,
                null,
                null);
    }

    @NotNull
    public static PurpleData load(@NotNull String tumorSample, @Nullable String referenceSample, @NotNull String qcFile,
            @NotNull String purityTsv, @NotNull String somaticDriverCatalogTsv, @NotNull String somaticVariantVcf,
            @NotNull String germlineDriverCatalogTsv, @NotNull String germlineVariantVcf, @Nullable String purpleGeneCopyNumberTsv,
            @Nullable String purpleSomaticCopynumberTsv, @Nullable RefGenomeVersion refGenomeVersion,
            @Nullable String alternativeTumorSampleId, @Nullable String alternativeReferenceSampleId) throws IOException {
        LOGGER.info("Loading PURPLE data from {}", new File(purityTsv).getParent());

        PurityContext purityContext = readPurityContext(qcFile, purityTsv);

        List<DriverCatalog> somaticDriverCatalog = DriverCatalogFile.read(somaticDriverCatalogTsv);
        LOGGER.info(" Loaded {} somatic driver catalog entries from {}", somaticDriverCatalog.size(), somaticDriverCatalogTsv);

        List<ReportableGainLoss> reportableGainsLosses = extractGainsLosses(somaticDriverCatalog);
        LOGGER.info("  Extracted {} reportable gains and losses from driver catalog", reportableGainsLosses.size());

        List<ReportableGainLoss> unreportedGainsLosses = Lists.newArrayList();
        if (purpleGeneCopyNumberTsv != null) {
            List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(purpleGeneCopyNumberTsv);
            LOGGER.info(" Loaded {} gene copy numbers entries from {}", geneCopyNumbers.size(), purpleGeneCopyNumberTsv);

            List<ReportableGainLoss> allGainsLosses =
                    extractAllGainsLosses(purityContext.qc().status(), purityContext.bestFit().ploidy(), geneCopyNumbers);
            LOGGER.debug("  Extracted {} gains and losses from gene copy numbers", allGainsLosses.size());

            unreportedGainsLosses = selectUnreportedGainsLosses(allGainsLosses, reportableGainsLosses);
            LOGGER.info("  Extracted {} additional unreported gains and losses", unreportedGainsLosses.size());
        }

        List<CnPerChromosomeArmData> cnPerChromosome = Lists.newArrayList();
        if (purpleSomaticCopynumberTsv != null && refGenomeVersion != null) {
            RefGenomeCoordinates refGenomeCoordinates =
                    refGenomeVersion == RefGenomeVersion.V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;
            cnPerChromosome = GenerateCnPerChromosome.fromPurpleSomaticCopynumberTsv(purpleSomaticCopynumberTsv, refGenomeCoordinates);
            LOGGER.info(" Loaded chromosomal arm copy numbers from {}", purpleSomaticCopynumberTsv);
        }

        List<ReportableVariant> reportableGermlineVariants = Lists.newArrayList();
        List<SomaticVariant> unreportedGermlineVariants = Lists.newArrayList();
        if (referenceSample != null || alternativeReferenceSampleId != null) {
            List<DriverCatalog> germlineDriverCatalog = DriverCatalogFile.read(germlineDriverCatalogTsv);
            LOGGER.info(" Loaded {} germline driver catalog entries from {}", germlineDriverCatalog.size(), germlineDriverCatalogTsv);

            List<SomaticVariant> germlineVariants = new SomaticVariantFactory().fromVCFFile(tumorSample,
                    referenceSample,
                    alternativeTumorSampleId,
                    alternativeReferenceSampleId,
                    germlineVariantVcf);
            reportableGermlineVariants = ReportableVariantFactory.toReportableGermlineVariants(germlineVariants, germlineDriverCatalog);
            LOGGER.info(" Loaded {} reportable germline variants from {}", reportableGermlineVariants.size(), germlineVariantVcf);

            unreportedGermlineVariants = selectUnreportedVariants(germlineVariants);
            LOGGER.info(" Loaded {} unreported germline variants from {}", unreportedGermlineVariants.size(), germlineVariantVcf);
        } else {
            LOGGER.info(" Skipped loading germline variants since no reference sample configured");
        }

        List<SomaticVariant> somaticVariants = SomaticVariantFactory.passOnlyInstance()
                .fromVCFFile(tumorSample, referenceSample, alternativeTumorSampleId, alternativeReferenceSampleId, somaticVariantVcf);
        List<ReportableVariant> reportableSomaticVariants =
                ReportableVariantFactory.toReportableSomaticVariants(somaticVariants, somaticDriverCatalog);
        LOGGER.info(" Loaded {} reportable somatic variants from {}", reportableSomaticVariants.size(), somaticVariantVcf);

        List<SomaticVariant> unreportedSomaticVariants = selectUnreportedVariants(somaticVariants);
        LOGGER.info(" Loaded {} unreported somatic variants from {}", unreportedSomaticVariants.size(), somaticVariantVcf);

        return ImmutablePurpleData.builder()
                .qc(purityContext.qc())
                .hasReliableQuality(purityContext.qc().pass())
                .fittedPurityMethod(purityContext.method())
                .wholeGenomeDuplication(purityContext.wholeGenomeDuplication())
                .hasReliablePurity(CheckPurpleQuality.checkHasReliablePurity(purityContext))
                .purity(purityContext.bestFit().purity())
                .minPurity(purityContext.score().minPurity())
                .maxPurity(purityContext.score().maxPurity())
                .ploidy(purityContext.bestFit().ploidy())
                .minPloidy(purityContext.score().minPloidy())
                .maxPloidy(purityContext.score().maxPloidy())
                .microsatelliteIndelsPerMb(purityContext.microsatelliteIndelsPerMb())
                .microsatelliteStatus(purityContext.microsatelliteStatus())
                .tumorMutationalBurdenPerMb(purityContext.tumorMutationalBurdenPerMb())
                .tumorMutationalLoad(purityContext.tumorMutationalLoad())
                .svTumorMutationalBurden(purityContext.svTumorMutationalBurden())
                .tumorMutationalLoadStatus(purityContext.tumorMutationalLoadStatus())
                .reportableSomaticVariants(reportableSomaticVariants)
                .unreportedSomaticVariants(unreportedSomaticVariants)
                .reportableGermlineVariants(reportableGermlineVariants)
                .unreportedGermlineVariants(unreportedGermlineVariants)
                .unreportedGainsLosses(unreportedGainsLosses)
                .reportableGainsLosses(reportableGainsLosses)
                .cnPerChromosome(cnPerChromosome)
                .build();
    }

    @NotNull
    private static PurityContext readPurityContext(@NotNull String qcFile, @NotNull String purityTsv) throws IOException {
        PurityContext purityContext = PurityContextFile.readWithQC(qcFile, purityTsv);

        DecimalFormat purityFormat = new DecimalFormat("#'%'");
        LOGGER.info("  QC status: {}", purityContext.qc().toString());
        LOGGER.info("  Tumor purity: {} ({}-{})",
                purityFormat.format(purityContext.bestFit().purity() * 100),
                purityFormat.format(purityContext.score().minPurity() * 100),
                purityFormat.format(purityContext.score().maxPurity() * 100));
        LOGGER.info("  Tumor ploidy: {} ({}-{})",
                purityContext.bestFit().ploidy(),
                purityContext.score().minPloidy(),
                purityContext.score().maxPloidy());
        LOGGER.info("  Fit method: {}", purityContext.method());
        LOGGER.info("  Whole genome duplication: {}", purityContext.wholeGenomeDuplication() ? "yes" : "no");
        LOGGER.info("  Microsatellite status: {}", purityContext.microsatelliteStatus().display());
        LOGGER.info("  Tumor mutational load status: {}", purityContext.tumorMutationalLoadStatus().display());

        return purityContext;
    }

    @NotNull
    private static List<SomaticVariant> selectUnreportedVariants(@NotNull List<SomaticVariant> variants) {
        List<SomaticVariant> filtered = Lists.newArrayList();
        for (SomaticVariant variant : variants) {
            if (!variant.reported()) {
                filtered.add(variant);
            }
        }
        return filtered;
    }

    @VisibleForTesting
    @NotNull
    static List<ReportableGainLoss> selectUnreportedGainsLosses(@NotNull List<ReportableGainLoss> allGainsLosses,
            @NotNull List<ReportableGainLoss> reportableGainsLosses) {
        List<ReportableGainLoss> unreportedGainsLosses = Lists.newArrayList();
        for (ReportableGainLoss gainLoss : allGainsLosses) {
            if (!reportableGainsLosses.contains(gainLoss)) {
                unreportedGainsLosses.add(gainLoss);
            }
        }
        return unreportedGainsLosses;
    }

    @NotNull
    private static List<ReportableGainLoss> extractGainsLosses(@NotNull List<DriverCatalog> drivers) {
        return drivers.stream()
                .filter(x -> x.driver() == DriverType.AMP || x.driver() == DriverType.PARTIAL_AMP || x.driver() == DriverType.DEL)
                .map(PurpleDataLoader::toReportableGainLoss)
                .collect(Collectors.toList());
    }

    @NotNull
    private static List<ReportableGainLoss> extractAllGainsLosses(@NotNull Set<PurpleQCStatus> qcStatus, double ploidy,
            @NotNull List<GeneCopyNumber> geneCopyNumbers) {
        List<DriverGene> allGenes = Lists.newArrayList();
        for (GeneCopyNumber geneCopyNumber : geneCopyNumbers) {
            allGenes.add(ImmutableDriverGene.builder()
                    .gene(geneCopyNumber.geneName())
                    .reportMissenseAndInframe(false)
                    .reportNonsenseAndFrameshift(false)
                    .reportSplice(false)
                    .reportDeletion(true)
                    .reportDisruption(false)
                    .reportAmplification(true)
                    .reportSomaticHotspot(false)
                    .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                    .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                    .reportGermlineDisruption(false)
                    .likelihoodType(DriverCategory.ONCO)
                    .build());
        }
        CNADrivers drivers = new CNADrivers(qcStatus, ImmutableDriverGenePanel.builder().driverGenes(allGenes).build());

        List<DriverCatalog> allGainLosses = Lists.newArrayList();
        allGainLosses.addAll(drivers.amplifications(ploidy, geneCopyNumbers));
        allGainLosses.addAll(drivers.deletions(geneCopyNumbers));

        return extractGainsLosses(allGainLosses);
    }

    @NotNull
    private static ReportableGainLoss toReportableGainLoss(@NotNull DriverCatalog driver) {
        return ImmutableReportableGainLoss.builder()
                .chromosome(driver.chromosome())
                .chromosomeBand(driver.chromosomeBand())
                .gene(driver.gene())
                .interpretation(CopyNumberInterpretation.fromCNADriver(driver))
                .minCopies(Math.round(Math.max(0, driver.minCopyNumber())))
                .maxCopies(Math.round(Math.max(0, driver.maxCopyNumber())))
                .build();
    }
}
