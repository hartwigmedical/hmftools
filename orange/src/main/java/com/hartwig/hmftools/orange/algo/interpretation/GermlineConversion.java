package com.hartwig.hmftools.orange.algo.interpretation;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.orange.algo.ImmutableOrangeReport;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.linx.ImmutableLinxInterpretedData;
import com.hartwig.hmftools.orange.algo.linx.LinxInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.ImmutablePurpleInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GermlineConversion {

    private static final Logger LOGGER = LogManager.getLogger(GermlineConversion.class);

    private GermlineConversion() {
    }

    @NotNull
    public static OrangeReport convertGermlineToSomatic(@NotNull OrangeReport report) {
        return ImmutableOrangeReport.builder()
                .from(report)
                .germlineMVLHPerGene(null)
                .purple(convertPurpleGermline(report.purple()))
                .linx(convertLinxGermline(report.linx()))
                .build();
    }

    @NotNull
    @VisibleForTesting
    static PurpleInterpretedData convertPurpleGermline(@NotNull PurpleInterpretedData purple) {
        // TODO Convert germline deletions into somatic deletions.
        return ImmutablePurpleInterpretedData.builder()
                .from(purple)
                .somaticDrivers(mergeGermlineDriversIntoSomatic(purple.somaticDrivers(), purple.germlineDrivers()))
                .germlineDrivers(null)
                .allSomaticVariants(mergeGermlineVariantsIntoSomatic(purple.allSomaticVariants(), purple.allGermlineVariants()))
                .reportableSomaticVariants(mergeGermlineVariantsIntoSomatic(purple.reportableSomaticVariants(),
                        purple.reportableGermlineVariants()))
                .additionalSuspectSomaticVariants(mergeGermlineVariantsIntoSomatic(purple.additionalSuspectSomaticVariants(),
                        purple.additionalSuspectGermlineVariants()))
                .allGermlineVariants(null)
                .reportableGermlineVariants(null)
                .additionalSuspectGermlineVariants(null)
                .allGermlineDeletions(null)
                .reportableGermlineDeletions(null)
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<DriverCatalog> mergeGermlineDriversIntoSomatic(@NotNull List<DriverCatalog> somaticDrivers,
            @NotNull List<DriverCatalog> germlineDrivers) {
        List<DriverCatalog> merged = Lists.newArrayList();
        for (DriverCatalog somaticDriver : somaticDrivers) {
            DriverCatalog matchingGermlineDriver = findMatchingGermlineDriver(somaticDriver, germlineDrivers);
            if (somaticDriver.driver() == DriverType.MUTATION && matchingGermlineDriver != null) {
                merged.add(mergeSomaticMutationWithGermline(somaticDriver, matchingGermlineDriver));
            } else {
                merged.add(somaticDriver);
            }
        }

        // TODO convert germline disruptions and germline deletions once their underlying data is converted.
        for (DriverCatalog germlineDriver : germlineDrivers) {
            if (germlineDriver.driver() == DriverType.GERMLINE_MUTATION
                    && findMatchingSomaticDriver(germlineDriver, somaticDrivers) == null) {
                merged.add(convertToSomaticDriver(germlineDriver));
            }
        }

        return merged;
    }

    @NotNull
    private static DriverCatalog mergeSomaticMutationWithGermline(@NotNull DriverCatalog somaticDriver,
            @NotNull DriverCatalog germlineDriver) {
        return ImmutableDriverCatalog.builder()
                .from(somaticDriver)
                .likelihoodMethod(LikelihoodMethod.HOTSPOT)
                .driverLikelihood(Math.max(somaticDriver.driverLikelihood(), germlineDriver.driverLikelihood()))
                .missense(somaticDriver.missense() + germlineDriver.missense())
                .nonsense(somaticDriver.nonsense() + germlineDriver.nonsense())
                .splice(somaticDriver.splice() + germlineDriver.splice())
                .frameshift(somaticDriver.frameshift() + germlineDriver.frameshift())
                .inframe(somaticDriver.inframe() + germlineDriver.inframe())
                .biallelic(somaticDriver.biallelic() || germlineDriver.biallelic())
                .build();
    }

    @Nullable
    private static DriverCatalog findMatchingGermlineDriver(@NotNull DriverCatalog somaticDriver,
            @NotNull List<DriverCatalog> germlineDrivers) {
        return find(germlineDrivers, DriverType.GERMLINE_MUTATION, somaticDriver.gene(), somaticDriver.transcript());
    }

    @Nullable
    private static DriverCatalog findMatchingSomaticDriver(@NotNull DriverCatalog germlineDriver,
            @NotNull List<DriverCatalog> somaticDrivers) {
        return find(somaticDrivers, DriverType.MUTATION, germlineDriver.gene(), germlineDriver.transcript());
    }

    @Nullable
    private static DriverCatalog find(@NotNull List<DriverCatalog> drivers, @NotNull DriverType driverTypeToFind,
            @NotNull String geneToFind, @NotNull String transcriptToFind) {
        for (DriverCatalog driver : drivers) {
            if (driver.driver() == driverTypeToFind && driver.gene().equals(geneToFind) && driver.transcript().equals(transcriptToFind)) {
                return driver;
            }
        }

        return null;
    }

    @NotNull
    private static DriverCatalog convertToSomaticDriver(@NotNull DriverCatalog germlineDriver) {
        if (!Doubles.equal(germlineDriver.driverLikelihood(), 1)) {
            LOGGER.warn("Germline driver converted to somatic with driver likelihood <> 1: {}", germlineDriver);
        }

        return ImmutableDriverCatalog.builder()
                .from(germlineDriver)
                .driver(DriverType.MUTATION)
                .likelihoodMethod(LikelihoodMethod.HOTSPOT)
                .build();
    }

    @NotNull
    private static List<PurpleVariant> mergeGermlineVariantsIntoSomatic(@NotNull List<PurpleVariant> somaticVariants,
            @NotNull List<PurpleVariant> germlineVariants) {
        List<PurpleVariant> merged = Lists.newArrayList();
        merged.addAll(somaticVariants);
        merged.addAll(germlineVariants);
        return merged;
    }

    @NotNull
    private static LinxInterpretedData convertLinxGermline(@NotNull LinxInterpretedData linx) {
        // TODO Convert germline disruptions to somatic disruptions.
        return ImmutableLinxInterpretedData.builder().from(linx).allGermlineDisruptions(null).reportableGermlineDisruptions(null).build();
    }
}
