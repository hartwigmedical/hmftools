package com.hartwig.hmftools.orange.algo.util;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.orange.algo.ImmutableOrangeReport;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.linx.ImmutableLinxInterpretedData;
import com.hartwig.hmftools.orange.algo.linx.LinxInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.ImmutablePurityPloidyFit;
import com.hartwig.hmftools.orange.algo.purple.ImmutablePurpleInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.PurityPloidyFit;
import com.hartwig.hmftools.orange.algo.purple.PurpleGainLoss;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GermlineConvertor {

    private static final Logger LOGGER = LogManager.getLogger(GermlineConvertor.class);

    @NotNull
    private final EnsemblDataCache ensemblDataCache;

    public GermlineConvertor(@NotNull final EnsemblDataCache ensemblDataCache) {
        this.ensemblDataCache = ensemblDataCache;
    }

    @NotNull
    public OrangeReport convertGermlineToSomatic(@NotNull OrangeReport report) {
        return ImmutableOrangeReport.builder()
                .from(report)
                .germlineMVLHPerGene(null)
                .purple(convertPurpleGermline(report.purple().fit().containsTumorCells(), report.purple(), report.linx()))
                .linx(convertLinxGermline(report.linx()))
                .build();
    }

    @NotNull
    @VisibleForTesting
    static PurpleInterpretedData convertPurpleGermline(boolean containsTumorCells, @NotNull PurpleInterpretedData purple,
            @NotNull LinxInterpretedData linx) {
        // In case tumor contains no tumor cells, we remove all germline events.
        List<DriverCatalog> mergedDrivers;
        List<PurpleVariant> additionalSomaticVariants;
        List<PurpleGainLoss> additionalSomaticGainsLosses;
        if (containsTumorCells) {
            mergedDrivers = mergeGermlineDriversIntoSomatic(purple.somaticDrivers(), purple.germlineDrivers());
            additionalSomaticVariants = toSomaticVariants(purple.reportableGermlineVariants());
            additionalSomaticGainsLosses = toSomaticGainLosses(purple.reportableGermlineDeletions());
        } else {
            mergedDrivers = purple.somaticDrivers();
            additionalSomaticVariants = Lists.newArrayList();
            additionalSomaticGainsLosses = Lists.newArrayList();
        }

        return ImmutablePurpleInterpretedData.builder()
                .from(purple)
                .fit(removeGermlineAberrations(purple.fit()))
                .somaticDrivers(mergedDrivers)
                .germlineDrivers(null)
                .addAllAllSomaticVariants(additionalSomaticVariants)
                .addAllReportableSomaticVariants(additionalSomaticVariants)
                .allGermlineVariants(null)
                .reportableGermlineVariants(null)
                .additionalSuspectGermlineVariants(null)
                .addAllAllSomaticGainsLosses(additionalSomaticGainsLosses)
                .addAllReportableSomaticGainsLosses(additionalSomaticGainsLosses)
                .allGermlineDeletions(null)
                .reportableGermlineDeletions(null)
                .build();
    }

    @NotNull
    private static PurityPloidyFit removeGermlineAberrations(@NotNull PurityPloidyFit fit) {
        return ImmutablePurityPloidyFit.builder()
                .from(fit)
                .qc(ImmutablePurpleQC.builder().from(fit.qc()).germlineAberrations(Sets.newHashSet()).build())
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<DriverCatalog> mergeGermlineDriversIntoSomatic(@NotNull List<DriverCatalog> somaticDrivers,
            @Nullable List<DriverCatalog> germlineDrivers) {
        List<DriverCatalog> merged = Lists.newArrayList();
        for (DriverCatalog somaticDriver : somaticDrivers) {
            DriverCatalog matchingGermlineDriver = findMatchingGermlineDriver(somaticDriver, germlineDrivers);
            if (somaticDriver.driver() == DriverType.MUTATION && matchingGermlineDriver != null) {
                merged.add(mergeSomaticMutationDriverWithGermline(somaticDriver, matchingGermlineDriver));
            } else {
                merged.add(somaticDriver);
            }
        }

        // TODO convert germline disruptions and germline deletions once their underlying data is converted.
        if (germlineDrivers != null) {
            for (DriverCatalog germlineDriver : germlineDrivers) {
                if (germlineDriver.driver() == DriverType.GERMLINE_MUTATION
                        && findMatchingSomaticDriver(germlineDriver, somaticDrivers) == null) {
                    merged.add(convertToSomaticDriver(germlineDriver));
                }
            }
        }

        return merged;
    }

    @NotNull
    private static List<PurpleVariant> toSomaticVariants(@Nullable List<PurpleVariant> reportableGermlineVariants) {
        return reportableGermlineVariants != null ? reportableGermlineVariants : Lists.newArrayList();
    }

    @NotNull
    private static List<PurpleGainLoss> toSomaticGainLosses(@Nullable List<GermlineDeletion> reportableGermlineDeletions) {
        if (reportableGermlineDeletions == null) {
            return Lists.newArrayList();
        }
        List<PurpleGainLoss> gainsLosses = Lists.newArrayList();
        for (GermlineDeletion deletion : reportableGermlineDeletions) {
            // TODO convert;
        }
        return gainsLosses;
    }

    @NotNull
    private static DriverCatalog mergeSomaticMutationDriverWithGermline(@NotNull DriverCatalog somaticDriver,
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
            @Nullable List<DriverCatalog> germlineDrivers) {
        if (germlineDrivers == null) {
            return null;
        }

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
    private static LinxInterpretedData convertLinxGermline(@NotNull LinxInterpretedData linx) {
        // TODO Convert germline disruptions to somatic disruptions.
        return ImmutableLinxInterpretedData.builder().from(linx).allGermlineDisruptions(null).reportableGermlineDisruptions(null).build();
    }
}
