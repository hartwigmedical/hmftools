package com.hartwig.hmftools.orange.algo.util;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.datamodel.linx.*;
import com.hartwig.hmftools.datamodel.purple.PurpleGainLoss;
import com.hartwig.hmftools.orange.algo.ImmutableOrangeReport;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.linx.ImmutableLinxInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.ImmutablePurityPloidyFit;
import com.hartwig.hmftools.orange.algo.purple.ImmutablePurpleInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.PurityPloidyFit;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariant;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;
import java.util.Map;

public final class GermlineConversion {

    private static final Logger LOGGER = LogManager.getLogger(GermlineConversion.class);

    private GermlineConversion() {
    }

    @NotNull
    public static OrangeReport convertGermlineToSomatic(@NotNull OrangeReport report) {
        boolean containsTumorCells = report.purple().fit().containsTumorCells();

        return ImmutableOrangeReport.builder()
                .from(report)
                .germlineMVLHPerGene(null)
                .purple(convertPurpleGermline(containsTumorCells, report.purple()))
                .linx(convertLinxGermline(containsTumorCells, report.linx()))
                .build();
    }

    @NotNull
    @VisibleForTesting
    static PurpleInterpretedData convertPurpleGermline(boolean containsTumorCells, @NotNull PurpleInterpretedData purple) {
        List<DriverCatalog> mergedDrivers;
        List<PurpleVariant> additionalReportableVariants;
        List<PurpleGainLoss> additionalReportableGainsLosses;
        if (containsTumorCells) {
            mergedDrivers = mergeGermlineDriversIntoSomatic(purple.somaticDrivers(), purple.germlineDrivers());
            additionalReportableVariants = toSomaticVariants(purple.reportableGermlineVariants());
            additionalReportableGainsLosses = toSomaticGainsLosses(purple.reportableGermlineFullLosses());
        } else {
            mergedDrivers = purple.somaticDrivers();
            additionalReportableVariants = Lists.newArrayList();
            additionalReportableGainsLosses = Lists.newArrayList();
        }

        return ImmutablePurpleInterpretedData.builder()
                .from(purple)
                .fit(removeGermlineAberrations(purple.fit()))
                .somaticDrivers(mergedDrivers)
                .germlineDrivers(null)
                .addAllAllSomaticVariants(additionalReportableVariants)
                .addAllReportableSomaticVariants(additionalReportableVariants)
                .allGermlineVariants(null)
                .reportableGermlineVariants(null)
                .additionalSuspectGermlineVariants(null)
                .addAllAllSomaticGainsLosses(additionalReportableGainsLosses)
                .addAllReportableSomaticGainsLosses(additionalReportableGainsLosses)
                .allGermlineDeletions(null)
                .allGermlineFullLosses(null)
                .reportableGermlineFullLosses(null)
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

        if (germlineDrivers != null) {
            for (DriverCatalog germlineDriver : germlineDrivers) {
                if (germlineDriver.driver() == DriverType.GERMLINE_MUTATION
                        && findMatchingSomaticDriver(germlineDriver, somaticDrivers, DriverType.MUTATION) == null) {
                    merged.add(convertToSomaticMutationDriver(germlineDriver));
                }

                if (germlineDriver.driver() == DriverType.GERMLINE_DELETION
                        && findMatchingSomaticDriver(germlineDriver, somaticDrivers, DriverType.DEL) == null) {
                    merged.add(convertToSomaticDeletionDriver(germlineDriver));
                }

                if (germlineDriver.driver() == DriverType.GERMLINE_DISRUPTION) {
                    merged.add(convertToSomaticDisruptionDriver(germlineDriver));
                }

                if (germlineDriver.driver() == DriverType.GERMLINE_HOM_DUP_DISRUPTION) {
                    merged.add(convertToSomaticHomozygousDisruptionDriver(germlineDriver));
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
    private static List<PurpleGainLoss> toSomaticGainsLosses(@Nullable List<PurpleGainLoss> reportableGermlineGainsLosses) {
        return reportableGermlineGainsLosses != null ? reportableGermlineGainsLosses : Lists.newArrayList();
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
            @NotNull List<DriverCatalog> somaticDrivers, @NotNull DriverType somaticTypeToFind) {
        return find(somaticDrivers, somaticTypeToFind, germlineDriver.gene(), germlineDriver.transcript());
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
    private static DriverCatalog convertToSomaticMutationDriver(@NotNull DriverCatalog germlineDriver) {
        if (!Doubles.equal(germlineDriver.driverLikelihood(), 1)) {
            LOGGER.warn("Germline driver converted to somatic with driver likelihood <> 1: {}", germlineDriver);
        }

        assert germlineDriver.driver() == DriverType.GERMLINE_MUTATION;

        return ImmutableDriverCatalog.builder()
                .from(germlineDriver)
                .driver(DriverType.MUTATION)
                .likelihoodMethod(LikelihoodMethod.HOTSPOT)
                .build();
    }

    @NotNull
    private static DriverCatalog convertToSomaticDeletionDriver(@NotNull DriverCatalog germlineDriver) {
        assert germlineDriver.driver() == DriverType.GERMLINE_DELETION;

        return ImmutableDriverCatalog.builder().from(germlineDriver).driver(DriverType.DEL).likelihoodMethod(LikelihoodMethod.DEL).build();
    }

    @NotNull
    private static DriverCatalog convertToSomaticDisruptionDriver(@NotNull DriverCatalog germlineDriver) {
        assert germlineDriver.driver() == DriverType.GERMLINE_DISRUPTION;

        return ImmutableDriverCatalog.builder()
                .from(germlineDriver)
                .driver(DriverType.DISRUPTION)
                .driverLikelihood(0D)
                .likelihoodMethod(LikelihoodMethod.DEL)
                .build();
    }

    @NotNull
    private static DriverCatalog convertToSomaticHomozygousDisruptionDriver(@NotNull DriverCatalog germlineDriver) {
        assert germlineDriver.driver() == DriverType.GERMLINE_HOM_DUP_DISRUPTION;

        return ImmutableDriverCatalog.builder()
                .from(germlineDriver)
                .driver(DriverType.HOM_DUP_DISRUPTION)
                .driverLikelihood(1D)
                .likelihoodMethod(LikelihoodMethod.DEL)
                .build();
    }

    @NotNull
    @VisibleForTesting
    static LinxRecord convertLinxGermline(boolean containsTumorCells, @NotNull LinxRecord linx) {
        List<LinxSvAnnotation> additionalStructuralVariants = Lists.newArrayList();
        List<LinxBreakend> additionalReportableBreakends = Lists.newArrayList();
        List<HomozygousDisruption> additionalHomozygousDisruptions = Lists.newArrayList();

        if (containsTumorCells) {
            Map<Integer, Integer> svIdMapping = buildSvIdMapping(linx.allSomaticStructuralVariants(), linx.allGermlineStructuralVariants());
            Map<Integer, Integer> clusterIdMapping =
                    buildClusterIdMapping(linx.allSomaticStructuralVariants(), linx.allGermlineStructuralVariants());
            Map<Integer, Integer> breakendIdMapping = buildBreakendIdMapping(linx.allSomaticBreakends(), linx.allGermlineBreakends());

            additionalStructuralVariants = toSomaticStructuralVariants(linx.allGermlineStructuralVariants(), svIdMapping, clusterIdMapping);
            additionalReportableBreakends = toSomaticBreakends(linx.reportableGermlineBreakends(), breakendIdMapping, svIdMapping);
            additionalHomozygousDisruptions = toSomaticHomozygousDisruptions(linx.germlineHomozygousDisruptions());
        }

        return ImmutableLinxRecord.builder()
                .from(linx)
                .addAllAllSomaticStructuralVariants(additionalStructuralVariants)
                .addAllAllSomaticBreakends(additionalReportableBreakends)
                .addAllReportableSomaticBreakends(additionalReportableBreakends)
                .addAllSomaticHomozygousDisruptions(additionalHomozygousDisruptions)
                .allGermlineStructuralVariants(null)
                .allGermlineBreakends(null)
                .reportableGermlineBreakends(null)
                .germlineHomozygousDisruptions(null)
                .build();
    }

    @NotNull
    private static List<LinxSvAnnotation> toSomaticStructuralVariants(@Nullable List<LinxSvAnnotation> germlineStructuralVariants,
            @NotNull Map<Integer, Integer> svIdMapping, @NotNull Map<Integer, Integer> clusterIdMapping) {
        List<LinxSvAnnotation> converted = Lists.newArrayList();
        if (germlineStructuralVariants != null) {
            for (LinxSvAnnotation structuralVariant : germlineStructuralVariants) {
                converted.add(ImmutableLinxSvAnnotation.builder()
                        .from(structuralVariant)
                        .svId(svIdMapping.get(structuralVariant.svId()))
                        .clusterId(clusterIdMapping.get(structuralVariant.clusterId()))
                        .build());
            }
        }
        return converted;
    }

    @NotNull
    private static List<LinxBreakend> toSomaticBreakends(@Nullable List<LinxBreakend> germlineBreakends,
            @NotNull Map<Integer, Integer> breakendIdMapping, @NotNull Map<Integer, Integer> svIdMapping) {
        List<LinxBreakend> converted = Lists.newArrayList();
        if (germlineBreakends != null) {
            for (LinxBreakend breakend : germlineBreakends) {
                converted.add(ImmutableLinxBreakend.builder()
                        .from(breakend)
                        .id(breakendIdMapping.get(breakend.id()))
                        .svId(svIdMapping.get(breakend.svId()))
                        .build());
            }
        }
        return converted;
    }

    @NotNull
    private static Map<Integer, Integer> buildSvIdMapping(@NotNull List<LinxSvAnnotation> allSomaticStructuralVariants,
            @Nullable List<LinxSvAnnotation> allGermlineStructuralVariants) {
        Map<Integer, Integer> svIdMapping = Maps.newHashMap();
        if (allGermlineStructuralVariants != null) {
            int newId = findMaxSvId(allSomaticStructuralVariants) + 1;
            for (LinxSvAnnotation structuralVariant : allGermlineStructuralVariants) {
                svIdMapping.put(structuralVariant.svId(), newId);
                newId++;
            }
        }
        return svIdMapping;
    }

    @VisibleForTesting
    static int findMaxSvId(@NotNull List<LinxSvAnnotation> structuralVariants) {
        int maxSvId = 0;
        for (LinxSvAnnotation structuralVariant : structuralVariants) {
            maxSvId = Math.max(maxSvId, structuralVariant.svId());
        }
        return maxSvId;
    }

    @NotNull
    private static Map<Integer, Integer> buildClusterIdMapping(@NotNull List<LinxSvAnnotation> allSomaticStructuralVariants,
            @Nullable List<LinxSvAnnotation> allGermlineStructuralVariants) {
        Map<Integer, Integer> clusterIdMapping = Maps.newHashMap();
        if (allGermlineStructuralVariants != null) {
            int newId = findMaxClusterId(allSomaticStructuralVariants) + 1;
            for (LinxSvAnnotation structuralVariant : allGermlineStructuralVariants) {
                clusterIdMapping.put(structuralVariant.clusterId(), newId);
                newId++;
            }
        }
        return clusterIdMapping;
    }

    @VisibleForTesting
    static int findMaxClusterId(@NotNull List<LinxSvAnnotation> structuralVariants) {
        int maxClusterId = 0;
        for (LinxSvAnnotation structuralVariant : structuralVariants) {
            maxClusterId = Math.max(maxClusterId, structuralVariant.clusterId());
        }
        return maxClusterId;
    }

    @NotNull
    private static Map<Integer, Integer> buildBreakendIdMapping(@NotNull List<LinxBreakend> allSomaticBreakends,
            @Nullable List<LinxBreakend> allGermlineBreakends) {
        Map<Integer, Integer> breakendIdMapping = Maps.newHashMap();
        if (allGermlineBreakends != null) {
            int newId = findMaxBreakendId(allSomaticBreakends) + 1;
            for (LinxBreakend breakend : allGermlineBreakends) {
                breakendIdMapping.put(breakend.id(), newId);
                newId++;
            }
        }
        return breakendIdMapping;
    }

    @VisibleForTesting
    static int findMaxBreakendId(@NotNull List<LinxBreakend> breakends) {
        int maxId = 0;
        for (LinxBreakend breakend : breakends) {
            maxId = Math.max(maxId, breakend.id());
        }
        return maxId;
    }

    @NotNull
    private static List<HomozygousDisruption> toSomaticHomozygousDisruptions(
            @Nullable List<HomozygousDisruption> germlineHomozygousDisruptions) {
        return germlineHomozygousDisruptions != null ? germlineHomozygousDisruptions : Lists.newArrayList();
    }
}
