package com.hartwig.hmftools.orange.algo.util;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.finding.GainDeletion;
import com.hartwig.hmftools.datamodel.finding.ImmutableGainDeletion;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleDriver;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleFit;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleVariant;
import com.hartwig.hmftools.datamodel.purple.ImmutableTumorStats;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleLikelihoodMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.TumorStats;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GermlineConversion
{
    @NotNull
    public static OrangeRecord convertGermlineToSomatic(final OrangeRecord report)
    {
        boolean containsTumorCells = report.purple().fit().containsTumorCells();

        PurpleRecord purple = convertPurpleGermline(containsTumorCells, report.purple());
        LinxRecord linx = convertLinxGermline(containsTumorCells, report.linx());

        return ImmutableOrangeRecord.builder()
                .from(report)
                .germlineMVLHPerGene(null)
                .purple(purple)
                .linx(linx)
                .build();
    }

    @NotNull
    @VisibleForTesting
    static PurpleRecord convertPurpleGermline(boolean containsTumorCells, final PurpleRecord purple)
    {
        TumorStats mergedTumorStats;
        List<PurpleDriver> mergedDrivers;
        List<PurpleVariant> additionalReportableVariants;
        List<GainDeletion> mergedAllSomaticGainsDels;
        List<GainDeletion> mergedReportableSomaticGainsDels;
        if(containsTumorCells)
        {
            mergedTumorStats = mergeTumorStats(purple);

            mergedDrivers = mergeGermlineDriversIntoSomatic(
                    purple.somaticDrivers(), purple.germlineDrivers(), purple.driverGermlineDeletions());

            additionalReportableVariants = toSomaticVariants(purple.driverGermlineVariants());

            mergedReportableSomaticGainsDels =
                    mergeGermlineFullDels(purple.driverGermlineDeletions(), purple.driverSomaticGainsDels());
        }
        else
        {
            mergedTumorStats = purple.tumorStats();
            mergedDrivers = purple.somaticDrivers();
            additionalReportableVariants = Lists.newArrayList();
            mergedReportableSomaticGainsDels = purple.driverSomaticGainsDels();
        }

        return ImmutablePurpleRecord.builder()
                .from(purple)
                .fit(removeGermlineAberrations(purple.fit()))
                .tumorStats(mergedTumorStats)
                .somaticDrivers(mergedDrivers)
                .germlineDrivers(null)
                .addAllOtherSomaticVariants(additionalReportableVariants)
                .addAllDriverSomaticVariants(additionalReportableVariants)
                .otherGermlineVariants(null)
                .driverGermlineVariants(null)
                .driverSomaticGainsDels(mergedReportableSomaticGainsDels)
                .otherGermlineDeletions(null)
                .driverGermlineDeletions(null)
                .allGermlineLossOfHeterozygosities(null)
                .driverGermlineLossOfHeterozygosities(null)
                .build();
    }

    static TumorStats mergeTumorStats(PurpleRecord purple)
    {
        return ImmutableTumorStats.builder()
                .from(purple.tumorStats())
                .hotspotMutationCount(
                        germlineHotspotCount(purple.driverGermlineVariants()) + purple.tumorStats().hotspotMutationCount())
                .build();
    }

    private static int germlineHotspotCount(@Nullable List<PurpleVariant> reportableGermlineVariants)
    {
        if(reportableGermlineVariants == null)
        {
            return 0;
        }

        return (int) reportableGermlineVariants.stream()
                .filter(variant -> variant.hotspot() == HotspotType.HOTSPOT)
                .count();
    }

    private static PurpleFit removeGermlineAberrations(final PurpleFit fit)
    {
        return ImmutablePurpleFit.builder()
                .from(fit)
                .qc(ImmutablePurpleQC.builder().from(fit.qc()).germlineAberrations(Sets.newHashSet()).build())
                .build();
    }

    @VisibleForTesting
    static List<PurpleDriver> mergeGermlineDriversIntoSomatic(final List<PurpleDriver> somaticDrivers,
            @Nullable List<PurpleDriver> germlineDrivers, @Nullable List<GainDeletion> reportableGermlineFullDels)
    {
        List<PurpleDriver> merged = Lists.newArrayList();
        for(PurpleDriver somaticDriver : somaticDrivers)
        {
            PurpleDriver matchingGermlineDriver = findMatchingGermlineDriver(somaticDriver, germlineDrivers);
            if(somaticDriver.type() == PurpleDriverType.MUTATION && matchingGermlineDriver != null)
            {
                merged.add(mergeSomaticMutationDriverWithGermline(somaticDriver, matchingGermlineDriver));
            }
            else
            {
                merged.add(somaticDriver);
            }
        }

        if(germlineDrivers != null)
        {
            for(PurpleDriver germlineDriver : germlineDrivers)
            {
                if(germlineDriver.type() == PurpleDriverType.GERMLINE_MUTATION
                        && findMatchingSomaticDriver(germlineDriver, somaticDrivers, PurpleDriverType.MUTATION) == null)
                {
                    merged.add(convertToSomaticMutationDriver(germlineDriver));
                }

                if(germlineDriver.type() == PurpleDriverType.GERMLINE_DELETION
                        && findMatchingSomaticDriver(germlineDriver, somaticDrivers, PurpleDriverType.DEL) == null
                        && reportableGermlineFullDels != null && reportableGermlineFullDels.stream()
                        .anyMatch(l -> l.gene().equals(germlineDriver.gene())))
                {
                    merged.add(convertToSomaticDeletionDriver(germlineDriver));
                }
            }
        }

        return merged;
    }

    @VisibleForTesting
    static List<PurpleVariant> toSomaticVariants(@Nullable List<PurpleVariant> reportableGermlineVariants)
    {
        return reportableGermlineVariants != null ? reportableGermlineVariants.stream()
                .map(variant -> variant.variantCopyNumber() >= 0.5
                        ? ImmutablePurpleVariant.builder().from(variant).subclonalLikelihood(0.0).build()
                        : ImmutablePurpleVariant.builder().from(variant).subclonalLikelihood(1.0).build())
                .collect(Collectors.toList()) : Lists.newArrayList();
    }

    private static List<GainDeletion> mergeGermlineFullDels(@Nullable List<GainDeletion> reportableGermlineFullDels,
            final List<GainDeletion> somaticGainDels)
    {
        if(reportableGermlineFullDels == null)
        {
            return somaticGainDels;
        }

        List<GainDeletion> mergedGainDels = Lists.newArrayList();
        for(GainDeletion somaticGainDel : somaticGainDels)
        {
            mergedGainDels.add(getGermlineCorrectedGainDel(somaticGainDel, reportableGermlineFullDels));
        }

        for(GainDeletion germlineDel : reportableGermlineFullDels)
        {
            boolean hasMatchingSomaticDel = somaticGainDels.stream().anyMatch(s -> delsMatch(germlineDel, s));
            if(!hasMatchingSomaticDel)
            {
                mergedGainDels.add(germlineDel);
            }
        }

        return mergedGainDels;
    }

    private static GainDeletion getGermlineCorrectedGainDel(final GainDeletion somaticGainDel,
            final List<GainDeletion> reportableGermlineFullDels)
    {
        if(!isDel(somaticGainDel))
        {
            return somaticGainDel;
        }

        List<GainDeletion> matchingGermlineDels =
                reportableGermlineFullDels.stream().filter(l -> delsMatch(l, somaticGainDel)).collect(Collectors.toList());
        if(matchingGermlineDels.isEmpty())
        {
            return somaticGainDel;
        }
        else if(matchingGermlineDels.size() == 1)
        {
            return mergeDels(somaticGainDel, matchingGermlineDels.get(0));
        }
        else
        {
            throw new IllegalStateException(
                    "More than one reportable germline del for transcript " + somaticGainDel.transcript() + " of gene "
                            + somaticGainDel.gene());
        }
    }

    private static boolean delsMatch(final GainDeletion germlineDel, final GainDeletion somaticGainDel)
    {
        return isDel(somaticGainDel) && germlineDel.gene().equals(somaticGainDel.gene()) && germlineDel.transcript()
                .equals(somaticGainDel.transcript());
    }

    private static GainDeletion mergeDels(final GainDeletion somaticDel, final GainDeletion germlineDel)
    {
        double minCopies = Math.min(somaticDel.minCopies(), germlineDel.minCopies());

        CopyNumberInterpretation interpretation;
        double maxCopies;
        if(somaticDel.interpretation() == germlineDel.interpretation())
        {
            interpretation = somaticDel.interpretation();
            maxCopies = Math.max(somaticDel.maxCopies(), germlineDel.maxCopies());
        }
        else if(isFullDel(somaticDel))
        {
            interpretation = CopyNumberInterpretation.FULL_DEL;
            maxCopies = somaticDel.maxCopies();
        }
        else
        {
            // germline full del and somatic partial del
            interpretation = CopyNumberInterpretation.FULL_DEL;
            maxCopies = germlineDel.maxCopies();
        }

        ReportedStatus reportedStatus = DriverUtils.maxReportedStatus(somaticDel.reportedStatus(), germlineDel.reportedStatus());
        DriverInterpretation driverInterpretation = DriverUtils.maxDriverInterpretation(somaticDel.driverInterpretation(), germlineDel.driverInterpretation());

        return ImmutableGainDeletion.builder()
                .findingKey(FindingKeys.gainDeletion(somaticDel.gene(), interpretation, somaticDel.isCanonical(), somaticDel.transcript()))
                .reportedStatus(reportedStatus)
                .driverInterpretation(driverInterpretation)
                .interpretation(interpretation)
                .chromosome(somaticDel.chromosome())
                .chromosomeBand(somaticDel.chromosomeBand())
                .gene(somaticDel.gene())
                .transcript(somaticDel.transcript())
                .isCanonical(somaticDel.isCanonical())
                .minCopies(minCopies)
                .maxCopies(maxCopies)
                .build();
    }

    private static boolean isDel(final GainDeletion gainDel)
    {
        return isFullDel(gainDel) || isPartialDel(gainDel);
    }

    private static boolean isFullDel(final GainDeletion gainDel)
    {
        return gainDel.interpretation() == CopyNumberInterpretation.FULL_DEL;
    }

    private static boolean isPartialDel(final GainDeletion gainDel)
    {
        return gainDel.interpretation() == CopyNumberInterpretation.PARTIAL_DEL;
    }

    private static PurpleDriver mergeSomaticMutationDriverWithGermline(final PurpleDriver somaticDriver,
            final PurpleDriver germlineDriver)
    {
        return ImmutablePurpleDriver.builder()
                .from(somaticDriver)
                .likelihoodMethod(PurpleLikelihoodMethod.HOTSPOT)
                .driverLikelihood(Math.max(somaticDriver.driverLikelihood(), germlineDriver.driverLikelihood()))
                .build();
    }

    @Nullable
    private static PurpleDriver findMatchingGermlineDriver(final PurpleDriver somaticDriver,
            @Nullable List<PurpleDriver> germlineDrivers)
    {
        if(germlineDrivers == null)
        {
            return null;
        }

        return find(germlineDrivers, PurpleDriverType.GERMLINE_MUTATION, somaticDriver.gene(), somaticDriver.transcript());
    }

    @Nullable
    private static PurpleDriver findMatchingSomaticDriver(final PurpleDriver germlineDriver, final List<PurpleDriver> somaticDrivers,
            final PurpleDriverType somaticTypeToFind)
    {
        return find(somaticDrivers, somaticTypeToFind, germlineDriver.gene(), germlineDriver.transcript());
    }

    @Nullable
    private static PurpleDriver find(final List<PurpleDriver> drivers, final PurpleDriverType PurpleDriverTypeToFind,
            final String geneToFind, final String transcriptToFind)
    {
        for(PurpleDriver driver : drivers)
        {
            if(driver.type() == PurpleDriverTypeToFind && driver.gene().equals(geneToFind) && driver.transcript()
                    .equals(transcriptToFind))
            {
                return driver;
            }
        }

        return null;
    }

    @NotNull
    private static PurpleDriver convertToSomaticMutationDriver(final PurpleDriver germlineDriver)
    {
        if(!Doubles.equal(germlineDriver.driverLikelihood(), 1))
        {
            LOGGER.warn("Germline driver converted to somatic with driver likelihood <> 1: {}", germlineDriver);
        }

        assert germlineDriver.type() == PurpleDriverType.GERMLINE_MUTATION;

        return ImmutablePurpleDriver.builder()
                .from(germlineDriver)
                .type(PurpleDriverType.MUTATION)
                .likelihoodMethod(PurpleLikelihoodMethod.HOTSPOT)
                .build();
    }

    @NotNull
    private static PurpleDriver convertToSomaticDeletionDriver(final PurpleDriver germlineDriver)
    {
        assert germlineDriver.type() == PurpleDriverType.GERMLINE_DELETION;

        return ImmutablePurpleDriver.builder()
                .from(germlineDriver)
                .type(PurpleDriverType.DEL)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .build();
    }

    @NotNull
    @VisibleForTesting
    static LinxRecord convertLinxGermline(boolean containsTumorCells, final LinxRecord linx)
    {
        List<LinxSvAnnotation> additionalStructuralVariants = Lists.newArrayList();
        List<LinxBreakend> additionalReportableBreakends = Lists.newArrayList();
        List<LinxHomozygousDisruption> additionalHomozygousDisruptions = Lists.newArrayList();

        if(containsTumorCells)
        {
            Map<Integer, Integer> svIdMapping = buildSvIdMapping(linx.allSomaticStructuralVariants(), linx.allGermlineStructuralVariants());
            Map<Integer, Integer> clusterIdMapping =
                    buildClusterIdMapping(linx.allSomaticStructuralVariants(), linx.allGermlineStructuralVariants());
            Map<Integer, Integer> breakendIdMapping = buildBreakendIdMapping(linx.otherSomaticBreakends(), linx.otherGermlineBreakends());

            additionalStructuralVariants = toSomaticStructuralVariants(linx.allGermlineStructuralVariants(), svIdMapping, clusterIdMapping);
            additionalReportableBreakends = toSomaticBreakends(linx.driverGermlineBreakends(), breakendIdMapping, svIdMapping);
            additionalHomozygousDisruptions = toSomaticHomozygousDisruptions(linx.germlineHomozygousDisruptions());
        }

        return ImmutableLinxRecord.builder()
                .from(linx)
                .addAllAllSomaticStructuralVariants(additionalStructuralVariants)
                .addAllOtherSomaticBreakends(additionalReportableBreakends)
                .addAllDriverSomaticBreakends(additionalReportableBreakends)
                .addAllSomaticHomozygousDisruptions(additionalHomozygousDisruptions)
                .allGermlineStructuralVariants(null)
                .otherGermlineBreakends(null)
                .driverGermlineBreakends(null)
                .germlineHomozygousDisruptions(null)
                .build();
    }

    @NotNull
    private static List<LinxSvAnnotation> toSomaticStructuralVariants(@Nullable List<LinxSvAnnotation> germlineStructuralVariants,
            final Map<Integer, Integer> svIdMapping, final Map<Integer, Integer> clusterIdMapping)
    {
        List<LinxSvAnnotation> converted = Lists.newArrayList();
        if(germlineStructuralVariants != null)
        {
            for(LinxSvAnnotation structuralVariant : germlineStructuralVariants)
            {
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
            final Map<Integer, Integer> breakendIdMapping, final Map<Integer, Integer> svIdMapping)
    {
        List<LinxBreakend> converted = Lists.newArrayList();
        if(germlineBreakends != null)
        {
            for(LinxBreakend breakend : germlineBreakends)
            {
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
    private static Map<Integer, Integer> buildSvIdMapping(final List<LinxSvAnnotation> allSomaticStructuralVariants,
            @Nullable List<LinxSvAnnotation> allGermlineStructuralVariants)
    {
        Map<Integer, Integer> svIdMapping = Maps.newHashMap();
        if(allGermlineStructuralVariants != null)
        {
            int newId = findMaxSvId(allSomaticStructuralVariants) + 1;
            for(LinxSvAnnotation structuralVariant : allGermlineStructuralVariants)
            {
                svIdMapping.put(structuralVariant.svId(), newId);
                newId++;
            }
        }
        return svIdMapping;
    }

    @VisibleForTesting
    static int findMaxSvId(final List<LinxSvAnnotation> structuralVariants)
    {
        int maxSvId = 0;
        for(LinxSvAnnotation structuralVariant : structuralVariants)
        {
            maxSvId = Math.max(maxSvId, structuralVariant.svId());
        }
        return maxSvId;
    }

    @NotNull
    private static Map<Integer, Integer> buildClusterIdMapping(final List<LinxSvAnnotation> allSomaticStructuralVariants,
            @Nullable List<LinxSvAnnotation> allGermlineStructuralVariants)
    {
        if(allGermlineStructuralVariants != null)
        {
            int idOffset = findMaxClusterId(allSomaticStructuralVariants);
            List<Integer> germlineClusterIds =
                    allGermlineStructuralVariants.stream().map(LinxSvAnnotation::clusterId).distinct().collect(Collectors.toList());
            return IntStream.range(0, germlineClusterIds.size())
                    .boxed()
                    .collect(Collectors.toMap(germlineClusterIds::get, clusterIndex -> clusterIndex + idOffset + 1));
        }
        else
        {
            return new HashMap<>();
        }
    }

    @VisibleForTesting
    static int findMaxClusterId(final List<LinxSvAnnotation> structuralVariants)
    {
        int maxClusterId = 0;
        for(LinxSvAnnotation structuralVariant : structuralVariants)
        {
            maxClusterId = Math.max(maxClusterId, structuralVariant.clusterId());
        }
        return maxClusterId;
    }

    @NotNull
    private static Map<Integer, Integer> buildBreakendIdMapping(final List<LinxBreakend> allSomaticBreakends,
            @Nullable List<LinxBreakend> allGermlineBreakends)
    {
        Map<Integer, Integer> breakendIdMapping = Maps.newHashMap();
        if(allGermlineBreakends != null)
        {
            int newId = findMaxBreakendId(allSomaticBreakends) + 1;
            for(LinxBreakend breakend : allGermlineBreakends)
            {
                breakendIdMapping.put(breakend.id(), newId);
                newId++;
            }
        }
        return breakendIdMapping;
    }

    @VisibleForTesting
    static int findMaxBreakendId(final List<LinxBreakend> breakends)
    {
        int maxId = 0;
        for(LinxBreakend breakend : breakends)
        {
            maxId = Math.max(maxId, breakend.id());
        }
        return maxId;
    }

    @NotNull
    private static List<LinxHomozygousDisruption> toSomaticHomozygousDisruptions(
            @Nullable List<LinxHomozygousDisruption> germlineHomozygousDisruptions)
    {
        return germlineHomozygousDisruptions != null ? germlineHomozygousDisruptions : Lists.newArrayList();
    }
}
