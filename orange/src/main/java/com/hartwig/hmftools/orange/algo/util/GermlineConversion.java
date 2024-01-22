package com.hartwig.hmftools.orange.algo.util;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.datamodel.linx.HomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleDriver;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleFit;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleHeterozygousDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleLikelihoodMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GermlineConversion {

    private static final Logger LOGGER = LogManager.getLogger(GermlineConversion.class);

    private GermlineConversion() {
    }

    @NotNull
    public static OrangeRecord convertGermlineToSomatic(@NotNull OrangeRecord report) {
        boolean containsTumorCells = report.purple().fit().containsTumorCells();

        return ImmutableOrangeRecord.builder()
                .from(report)
                .germlineMVLHPerGene(null)
                .purple(convertPurpleGermline(containsTumorCells, report.purple()))
                .linx(convertLinxGermline(containsTumorCells, report.linx()))
                .build();
    }

    @NotNull
    @VisibleForTesting
    static PurpleRecord convertPurpleGermline(boolean containsTumorCells, @NotNull PurpleRecord purple)
    {
        List<PurpleDriver> mergedDrivers;
        List<PurpleVariant> additionalReportableVariants;
        List<PurpleGainLoss> additionalReportableGainsLosses;
        List<PurpleGeneCopyNumber> additionalSuspectGeneCopyNumbersWithLOH;
        if(containsTumorCells)
        {
            mergedDrivers = mergeGermlineDriversIntoSomatic(purple.somaticDrivers(), purple.germlineDrivers());
            additionalReportableVariants = toSomaticVariants(purple.reportableGermlineVariants());
            additionalReportableGainsLosses = toSomaticGainsLosses(purple.reportableGermlineFullLosses());
            additionalSuspectGeneCopyNumbersWithLOH = getAdditionalSuspectGeneCopyNumberWithLOH(purple);
        }
        else
        {
            mergedDrivers = purple.somaticDrivers();
            additionalReportableVariants = Lists.newArrayList();
            additionalReportableGainsLosses = Lists.newArrayList();
            additionalSuspectGeneCopyNumbersWithLOH = Lists.newArrayList();
        }

        return ImmutablePurpleRecord.builder()
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
                .addAllSuspectGeneCopyNumbersWithLOH(additionalSuspectGeneCopyNumbersWithLOH)
                .allGermlineHeterozygousDeletions(null)
                .reportableGermlineHeterozygousDeletions(null)
                .build();
    }

    @NotNull
    private static PurpleFit removeGermlineAberrations(@NotNull PurpleFit fit) {
        return ImmutablePurpleFit.builder()
                .from(fit)
                .qc(ImmutablePurpleQC.builder().from(fit.qc()).germlineAberrations(Sets.newHashSet()).build())
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<PurpleDriver> mergeGermlineDriversIntoSomatic(@NotNull List<PurpleDriver> somaticDrivers,
            @Nullable List<PurpleDriver> germlineDrivers) {
        List<PurpleDriver> merged = Lists.newArrayList();
        for (PurpleDriver somaticDriver : somaticDrivers) {
            PurpleDriver matchingGermlineDriver = findMatchingGermlineDriver(somaticDriver, germlineDrivers);
            if (somaticDriver.driver() == PurpleDriverType.MUTATION && matchingGermlineDriver != null) {
                merged.add(mergeSomaticMutationDriverWithGermline(somaticDriver, matchingGermlineDriver));
            } else {
                merged.add(somaticDriver);
            }
        }

        if (germlineDrivers != null) {
            for (PurpleDriver germlineDriver : germlineDrivers) {
                if (germlineDriver.driver() == PurpleDriverType.GERMLINE_MUTATION
                        && findMatchingSomaticDriver(germlineDriver, somaticDrivers, PurpleDriverType.MUTATION) == null) {
                    merged.add(convertToSomaticMutationDriver(germlineDriver));
                }

                if (germlineDriver.driver() == PurpleDriverType.GERMLINE_DELETION
                        && findMatchingSomaticDriver(germlineDriver, somaticDrivers, PurpleDriverType.DEL) == null) {
                    merged.add(convertToSomaticDeletionDriver(germlineDriver));
                }

                if (germlineDriver.driver() == PurpleDriverType.GERMLINE_DISRUPTION) {
                    merged.add(convertToSomaticDisruptionDriver(germlineDriver));
                }

                if (germlineDriver.driver() == PurpleDriverType.GERMLINE_HOM_DUP_DISRUPTION) {
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
    private static List<PurpleGeneCopyNumber> getAdditionalSuspectGeneCopyNumberWithLOH(@NotNull PurpleRecord purple)
    {
        List<PurpleHeterozygousDeletion> reportableGermlineHeterozygousDeletions = purple.reportableGermlineHeterozygousDeletions();
        if(reportableGermlineHeterozygousDeletions == null)
        {
            return Lists.newArrayList();
        }

        List<String> relevantGenes = reportableGermlineHeterozygousDeletions.stream()
                .map(PurpleHeterozygousDeletion::gene)
                .distinct()
                .collect(Collectors.toList());

        List<PurpleGeneCopyNumber> additionalSuspectGeneCopyNumberWithLOH = Lists.newArrayList();
        for(String gene : relevantGenes)
        {
            PurpleGeneCopyNumber candidate = findMatchingGeneCopyNumber(gene, purple.allSomaticGeneCopyNumbers());

            if(candidate != null
                    && shouldAddSuspectGeneCopyNumberWithLOH(candidate, purple.suspectGeneCopyNumbersWithLOH(), purple.reportableGermlineFullLosses()))
            {
                PurpleGeneCopyNumber adjustedGeneCopyNumber = correctForGermlineImpact(candidate, reportableGermlineHeterozygousDeletions);
                additionalSuspectGeneCopyNumberWithLOH.add(adjustedGeneCopyNumber);
            }
        }
        return additionalSuspectGeneCopyNumberWithLOH;
    }

    private static boolean shouldAddSuspectGeneCopyNumberWithLOH(@NotNull PurpleGeneCopyNumber candidateSuspectGeneCopyNumber,
            @NotNull List<PurpleGeneCopyNumber> suspectGeneCopyNumbersWithLOH, @Nullable List<PurpleGainLoss> reportableGermlineFullLosses)
    {
        String gene = candidateSuspectGeneCopyNumber.geneName();

        boolean alreadySuspectGeneCopyNumberWithLOH = findMatchingGeneCopyNumber(gene, suspectGeneCopyNumbersWithLOH) != null;
        boolean hasReportableGermlineFullLoss =
                reportableGermlineFullLosses != null && reportableGermlineFullLosses.stream().anyMatch(l -> l.gene().equals(gene));
        boolean notFullyLostInTumor = notFullyLostInTumor(candidateSuspectGeneCopyNumber);

        return !alreadySuspectGeneCopyNumberWithLOH && !hasReportableGermlineFullLoss && notFullyLostInTumor;
    }

    @Nullable
    private static PurpleGeneCopyNumber findMatchingGeneCopyNumber(@NotNull String gene,
            @NotNull List<PurpleGeneCopyNumber> geneCopyNumbers)
    {
        for(PurpleGeneCopyNumber geneCopyNumber : geneCopyNumbers)
        {
            if(geneCopyNumber.geneName().equals(gene))
            {
                return geneCopyNumber;
            }
        }
        return null;
    }

    private static boolean notFullyLostInTumor(@NotNull PurpleGeneCopyNumber geneCopyNumber)
    {
        return geneCopyNumber.minCopyNumber() > 0.5;
    }

    @NotNull
    private static PurpleGeneCopyNumber correctForGermlineImpact(@NotNull PurpleGeneCopyNumber geneCopyNumber,
            @NotNull List<PurpleHeterozygousDeletion> reportableGermlineHeterozygousDeletions)
    {
        double adjustedMinCopyNumber = adjustMinCopyNumberForGermlineImpact(geneCopyNumber, reportableGermlineHeterozygousDeletions);
        return ImmutablePurpleGeneCopyNumber.builder()
                .from(geneCopyNumber)
                .minCopyNumber(adjustedMinCopyNumber)
                .minMinorAlleleCopyNumber(0D)
                .build();
    }

    private static double adjustMinCopyNumberForGermlineImpact(@NotNull PurpleGeneCopyNumber geneCopyNumber,
            @NotNull List<PurpleHeterozygousDeletion> reportableGermlineHeterozygousDeletions)
    {
        OptionalDouble minimumTumorCopyNumber = reportableGermlineHeterozygousDeletions.stream()
                .filter(d -> d.gene().equals(geneCopyNumber.geneName()))
                .mapToDouble(PurpleHeterozygousDeletion::minCopies)
                .min();
        if(minimumTumorCopyNumber.isPresent())
        {
            return Math.min(minimumTumorCopyNumber.getAsDouble(), geneCopyNumber.minCopyNumber());
        }
        else
        {
            return geneCopyNumber.minCopyNumber();
        }
    }

    @NotNull
    private static PurpleDriver mergeSomaticMutationDriverWithGermline(@NotNull PurpleDriver somaticDriver,
            @NotNull PurpleDriver germlineDriver) {
        return ImmutablePurpleDriver.builder()
                .from(somaticDriver)
                .likelihoodMethod(PurpleLikelihoodMethod.HOTSPOT)
                .driverLikelihood(Math.max(somaticDriver.driverLikelihood(), germlineDriver.driverLikelihood()))
                .build();
    }

    @Nullable
    private static PurpleDriver findMatchingGermlineDriver(@NotNull PurpleDriver somaticDriver,
            @Nullable List<PurpleDriver> germlineDrivers) {
        if (germlineDrivers == null) {
            return null;
        }

        return find(germlineDrivers, PurpleDriverType.GERMLINE_MUTATION, somaticDriver.gene(), somaticDriver.transcript());
    }

    @Nullable
    private static PurpleDriver findMatchingSomaticDriver(@NotNull PurpleDriver germlineDriver, @NotNull List<PurpleDriver> somaticDrivers,
            @NotNull PurpleDriverType somaticTypeToFind) {
        return find(somaticDrivers, somaticTypeToFind, germlineDriver.gene(), germlineDriver.transcript());
    }

    @Nullable
    private static PurpleDriver find(@NotNull List<PurpleDriver> drivers, @NotNull PurpleDriverType PurpleDriverTypeToFind,
            @NotNull String geneToFind, @NotNull String transcriptToFind) {
        for (PurpleDriver driver : drivers) {
            if (driver.driver() == PurpleDriverTypeToFind && driver.gene().equals(geneToFind) && driver.transcript()
                    .equals(transcriptToFind)) {
                return driver;
            }
        }

        return null;
    }

    @NotNull
    private static PurpleDriver convertToSomaticMutationDriver(@NotNull PurpleDriver germlineDriver) {
        if (!Doubles.equal(germlineDriver.driverLikelihood(), 1)) {
            LOGGER.warn("Germline driver converted to somatic with driver likelihood <> 1: {}", germlineDriver);
        }

        assert germlineDriver.driver() == PurpleDriverType.GERMLINE_MUTATION;

        return ImmutablePurpleDriver.builder()
                .from(germlineDriver)
                .driver(PurpleDriverType.MUTATION)
                .likelihoodMethod(PurpleLikelihoodMethod.HOTSPOT)
                .build();
    }

    @NotNull
    private static PurpleDriver convertToSomaticDeletionDriver(@NotNull PurpleDriver germlineDriver) {
        assert germlineDriver.driver() == PurpleDriverType.GERMLINE_DELETION;

        return ImmutablePurpleDriver.builder()
                .from(germlineDriver)
                .driver(PurpleDriverType.DEL)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .build();
    }

    @NotNull
    private static PurpleDriver convertToSomaticDisruptionDriver(@NotNull PurpleDriver germlineDriver) {
        assert germlineDriver.driver() == PurpleDriverType.GERMLINE_DISRUPTION;

        return ImmutablePurpleDriver.builder()
                .from(germlineDriver)
                .driver(PurpleDriverType.DISRUPTION)
                .driverLikelihood(0D)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .build();
    }

    @NotNull
    private static PurpleDriver convertToSomaticHomozygousDisruptionDriver(@NotNull PurpleDriver germlineDriver) {
        assert germlineDriver.driver() == PurpleDriverType.GERMLINE_HOM_DUP_DISRUPTION;

        return ImmutablePurpleDriver.builder()
                .from(germlineDriver)
                .driver(PurpleDriverType.HOM_DUP_DISRUPTION)
                .driverLikelihood(1D)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
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
        if (allGermlineStructuralVariants != null) {
            int idOffset = findMaxClusterId(allSomaticStructuralVariants);
            List<Integer> germlineClusterIds =
                    allGermlineStructuralVariants.stream().map(LinxSvAnnotation::clusterId).distinct().collect(Collectors.toList());
            return IntStream.range(0, germlineClusterIds.size())
                    .boxed()
                    .collect(Collectors.toMap(germlineClusterIds::get, clusterIndex -> clusterIndex + idOffset + 1));
        } else {
            return new HashMap<>();
        }
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
