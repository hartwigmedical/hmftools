package com.hartwig.hmftools.protect.purple;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ImmutableReportableGainLoss;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.protect.ProtectConfig;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class PurpleDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(PurpleDataLoader.class);

    private PurpleDataLoader() {
    }

    @NotNull
    public static PurpleData load(@NotNull ProtectConfig config) throws IOException {
        return load(config.tumorSampleId(),
                config.purpleQcFile(),
                config.purplePurityTsv(),
                config.purpleDriverCatalogTsv(),
                config.purpleSomaticVariantVcf());
    }

    @NotNull
    private static PurpleData load(@NotNull String sample, @NotNull String qcFile, @NotNull String purityFile,
            @NotNull String driverCatalogFile, @NotNull String somaticVcf) throws IOException {
        LOGGER.info("Loading PURPLE data from {}", new File(purityFile).getParent());

        PurityContext purityContext = PurityContextFile.readWithQC(qcFile, purityFile);
        LOGGER.info("  QC status: {}", purityContext.qc().toString());
        LOGGER.info("  Tumor purity: {}", new DecimalFormat("#'%'").format(purityContext.bestFit().purity() * 100));
        LOGGER.info("  Fit method: {}", purityContext.method());
        LOGGER.info("  Average tumor ploidy: {}", purityContext.bestFit().ploidy());
        LOGGER.info("  Whole genome duplication: {}", purityContext.wholeGenomeDuplication() ? "yes" : "no");
        LOGGER.info("  Microsatellite status: {}", purityContext.microsatelliteStatus().display());
        LOGGER.info("  Tumor mutational load status: {}", purityContext.tumorMutationalLoadStatus().display());

        List<DriverCatalog> driverCatalog = DriverCatalogFile.read(driverCatalogFile);
        LOGGER.info(" Loaded {} driver catalog entries from {}", driverCatalog.size(), driverCatalogFile);

        List<ReportableGainLoss> copyNumberAlterations = driverCatalog.stream()
                .filter(x -> x.driver() == DriverType.AMP || x.driver() == DriverType.DEL)
                .map(PurpleDataLoader::copyNumberAlteration)
                .collect(Collectors.toList());
        LOGGER.info("  Extracted {} reportable copy number alterations from driver catalog", copyNumberAlterations.size());

        List<SomaticVariant> variants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(sample, somaticVcf);
        List<ReportableVariant> reportableVariants = ReportableVariantFactory.reportableSomaticVariants(variants, driverCatalog);
        LOGGER.info(" Loaded {} reportable somatic variants from {}", reportableVariants.size(), somaticVcf);

        return ImmutablePurpleData.builder()
                .purity(purityContext.bestFit().purity())
                .ploidy(purityContext.bestFit().ploidy())
                .microsatelliteIndelsPerMb(purityContext.microsatelliteIndelsPerMb())
                .microsatelliteStatus(purityContext.microsatelliteStatus())
                .tumorMutationalBurdenPerMb(purityContext.tumorMutationalBurdenPerMb())
                .tumorMutationalLoad(purityContext.tumorMutationalLoad())
                .tumorMutationalLoadStatus(purityContext.tumorMutationalLoadStatus())
                .somaticVariants(reportableVariants)
                .copyNumberAlterations(copyNumberAlterations)
                .build();
    }

    @NotNull
    private static ReportableGainLoss copyNumberAlteration(@NotNull DriverCatalog driver) {
        return ImmutableReportableGainLoss.builder()
                .chromosome(driver.chromosome())
                .chromosomeBand(driver.chromosomeBand())
                .gene(driver.gene())
                .interpretation(CopyNumberInterpretation.fromCNADriver(driver))
                .copies(Math.round(Math.max(0, driver.minCopyNumber())))
                .build();
    }
}
