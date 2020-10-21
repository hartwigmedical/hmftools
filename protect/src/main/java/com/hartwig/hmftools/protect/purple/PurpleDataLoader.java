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
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.protect.ProtectConfig;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PurpleDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(PurpleDataLoader.class);

    private PurpleDataLoader() {
    }

    @NotNull
    public static PurpleData load(@NotNull ProtectConfig config) throws IOException {
        return load(config.tumorSampleId(),
                config.purpleQcFile(),
                config.purplePurityTsv(),
                config.purpleGeneCnvTsv(),
                config.purpleDriverCatalogTsv(),
                config.purpleSomaticVariantVcf());
    }

    @NotNull
    public static PurpleData load(String sample, String qcFile, String purityFile, String geneCopyNumberFile, String driverCatalogFile,
            String somaticVcf) throws IOException {
        final PurityContext purityContext = PurityContextFile.readWithQC(qcFile, purityFile);
        LOGGER.info("Loaded PURPLE data from {}", new File(purityFile).getParent());
        LOGGER.info(" Purity: {}", new DecimalFormat("#'%'").format(purityContext.bestFit().purity() * 100));
        LOGGER.info(" Average tumor ploidy: {}", purityContext.bestFit().ploidy());
        LOGGER.info(" Fit method: {}", purityContext.method());
        LOGGER.info(" Whole genome duplication: {}", purityContext.wholeGenomeDuplication() ? "yes" : "no");

        PurpleQC purpleQC = purityContext.qc();
        LOGGER.info(" QC status: {}", purpleQC.toString());

        //TODO: Hopefully get rid of this
        List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(geneCopyNumberFile);
        LOGGER.info(" Gene copy numbers: {}", geneCopyNumbers.size());

        //TODO: Move remainder into purple
        List<DriverCatalog> driverCatalog = DriverCatalogFile.read(driverCatalogFile);
        List<SomaticVariant> variants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(sample, somaticVcf);
        List<ReportableVariant> reportableVariants = ReportableVariantFactory.reportableSomaticVariants(driverCatalog, variants);
        LOGGER.info(" Reportable somatic variants: {}", reportableVariants.size());

        List<ReportableGainLoss> copyNumberAlterations = driverCatalog.stream()
                .filter(x -> x.driver() == DriverType.AMP || x.driver() == DriverType.DEL)
                .map(PurpleDataLoader::copyNumberAlteration)
                .collect(Collectors.toList());
        LOGGER.info(" Reportable copy number alterations: {}", copyNumberAlterations.size());

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
                .geneCopyNumbers(geneCopyNumbers)
                .build();

    }

    @NotNull
    private static ReportableGainLoss copyNumberAlteration(DriverCatalog driver) {
        return ImmutableReportableGainLoss.builder()
                .chromosome(driver.chromosome())
                .chromosomeBand(driver.chromosomeBand())
                .gene(driver.gene())
                .interpretation(CopyNumberInterpretation.fromCNADriver(driver))
                .copies(Math.round(Math.max(0, driver.minCopyNumber())))
                .build();
    }
}
