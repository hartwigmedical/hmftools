package com.hartwig.hmftools.orange.algo.util;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.Doubles;
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
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleVariant;
import com.hartwig.hmftools.datamodel.purple.ImmutableTumorStats;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleFittedPurityMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleLikelihoodMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleLossOfHeterozygosity;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.TumorStats;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GermlineConversion
{
    @NotNull
    public static OrangeRecord convertGermlineToSomatic(@NotNull OrangeRecord report)
    {
        boolean containsTumorCells = containsTumorCells(report.purple().fit());

        return ImmutableOrangeRecord.builder()
                .from(report)
                .germlineMVLHPerGene(null)
                .purple(convertPurpleGermline(containsTumorCells, report.purple()))
                .linx(convertLinxGermline(containsTumorCells, report.linx()))
                .build();
    }

    private static boolean containsTumorCells(@NotNull PurpleFit purpleFit)
    {
        return purpleFit.fittedPurityMethod() != PurpleFittedPurityMethod.NO_TUMOR
                && !purpleFit.qc().status().contains(PurpleQCStatus.FAIL_NO_TUMOR);
    }

    @NotNull
    @VisibleForTesting
    static PurpleRecord convertPurpleGermline(boolean containsTumorCells, @NotNull PurpleRecord purple)
    {
        TumorStats mergedTumorStats;
        List<PurpleDriver> mergedDrivers;
        List<PurpleVariant> additionalReportableVariants;
        List<PurpleGainDeletion> mergedAllSomaticGainsDels;
        List<PurpleGainDeletion> mergedReportableSomaticGainsDels;
        List<PurpleGeneCopyNumber> mergedSuspectGeneCopyNumbersWithLOH;
        if(containsTumorCells)
        {
            mergedTumorStats = mergeTumorStats(purple);
            mergedDrivers =
                    mergeGermlineDriversIntoSomatic(purple.somaticDrivers(), purple.germlineDrivers(),
                            purple.reportableGermlineFullDels());
            additionalReportableVariants = toSomaticVariants(purple.reportableGermlineVariants());
            mergedAllSomaticGainsDels = mergeGermlineFullDels(purple.reportableGermlineFullDels(), purple.allSomaticGainsDels());
            mergedReportableSomaticGainsDels =
                    mergeGermlineFullDels(purple.reportableGermlineFullDels(), purple.reportableSomaticGainsDels());
            mergedSuspectGeneCopyNumbersWithLOH = mergeSuspectGeneCopyNumberWithLOH(purple);
        }
        else
        {
            mergedTumorStats = purple.tumorStats();
            mergedDrivers = purple.somaticDrivers();
            additionalReportableVariants = Lists.newArrayList();
            mergedAllSomaticGainsDels = purple.allSomaticGainsDels();
            mergedReportableSomaticGainsDels = purple.reportableSomaticGainsDels();
            mergedSuspectGeneCopyNumbersWithLOH = purple.suspectGeneCopyNumbersWithLOH();
        }

        return ImmutablePurpleRecord.builder()
                .from(purple)
                .fit(removeGermlineAberrations(purple.fit()))
                .tumorStats(mergedTumorStats)
                .somaticDrivers(mergedDrivers)
                .germlineDrivers(null)
                .addAllAllSomaticVariants(additionalReportableVariants)
                .addAllReportableSomaticVariants(additionalReportableVariants)
                .allGermlineVariants(null)
                .reportableGermlineVariants(null)
                .additionalSuspectGermlineVariants(null)
                .allSomaticGainsDels(mergedAllSomaticGainsDels)
                .reportableSomaticGainsDels(mergedReportableSomaticGainsDels)
                .allGermlineDeletions(null)
                .allGermlineFullDels(null)
                .reportableGermlineFullDels(null)
                .suspectGeneCopyNumbersWithLOH(mergedSuspectGeneCopyNumbersWithLOH)
                .allGermlineLossOfHeterozygosities(null)
                .reportableGermlineLossOfHeterozygosities(null)
                .build();
    }

    @NotNull
    static TumorStats mergeTumorStats(PurpleRecord purple)
    {
        return ImmutableTumorStats.builder()
                .from(purple.tumorStats())
                .hotspotMutationCount(
                        germlineHotspotCount(purple.reportableGermlineVariants()) + purple.tumorStats().hotspotMutationCount())
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

    @NotNull
    private static PurpleFit removeGermlineAberrations(@NotNull PurpleFit fit)
    {
        return ImmutablePurpleFit.builder()
                .from(fit)
                .qc(ImmutablePurpleQC.builder().from(fit.qc()).germlineAberrations(Sets.newHashSet()).build())
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<PurpleDriver> mergeGermlineDriversIntoSomatic(@NotNull List<PurpleDriver> somaticDrivers,
            @Nullable List<PurpleDriver> germlineDrivers, @Nullable List<PurpleGainDeletion> reportableGermlineFullDels)
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

    @NotNull
    @VisibleForTesting
    static List<PurpleVariant> toSomaticVariants(@Nullable List<PurpleVariant> reportableGermlineVariants)
    {
        return reportableGermlineVariants != null ? reportableGermlineVariants.stream()
                .map(variant -> variant.variantCopyNumber() >= 0.5
                        ? ImmutablePurpleVariant.builder().from(variant).subclonalLikelihood(0.0).build()
                        : ImmutablePurpleVariant.builder().from(variant).subclonalLikelihood(1.0).build())
                .collect(Collectors.toList()) : Lists.newArrayList();
    }

    @NotNull
    private static List<PurpleGainDeletion> mergeGermlineFullDels(@Nullable List<PurpleGainDeletion> reportableGermlineFullDels,
            @NotNull List<PurpleGainDeletion> somaticGainDels)
    {
        if(reportableGermlineFullDels == null)
        {
            return somaticGainDels;
        }

        List<PurpleGainDeletion> mergedGainDels = Lists.newArrayList();
        for(PurpleGainDeletion somaticGainDel : somaticGainDels)
        {
            mergedGainDels.add(getGermlineCorrectedGainDel(somaticGainDel, reportableGermlineFullDels));
        }

        for(PurpleGainDeletion germlineDel : reportableGermlineFullDels)
        {
            boolean hasMatchingSomaticDel = somaticGainDels.stream().anyMatch(s -> DelsMatch(germlineDel, s));
            if(!hasMatchingSomaticDel)
            {
                mergedGainDels.add(germlineDel);
            }
        }

        return mergedGainDels;
    }

    @NotNull
    private static PurpleGainDeletion getGermlineCorrectedGainDel(@NotNull PurpleGainDeletion somaticGainDel,
            @NotNull List<PurpleGainDeletion> reportableGermlineFullDels)
    {
        if(!isDel(somaticGainDel))
        {
            return somaticGainDel;
        }

        List<PurpleGainDeletion> matchingGermlineDels =
                reportableGermlineFullDels.stream().filter(l -> DelsMatch(l, somaticGainDel)).collect(Collectors.toList());
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

    private static boolean DelsMatch(@NotNull PurpleGainDeletion germlineDel, @NotNull PurpleGainDeletion somaticGainDel)
    {
        return isDel(somaticGainDel) && germlineDel.gene().equals(somaticGainDel.gene()) && germlineDel.transcript()
                .equals(somaticGainDel.transcript());
    }

    @NotNull
    private static PurpleGainDeletion mergeDels(@NotNull PurpleGainDeletion somaticDel, @NotNull PurpleGainDeletion germlineDel)
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

        return ImmutablePurpleGainDeletion.builder()
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

    private static boolean isDel(@NotNull PurpleGainDeletion gainDel)
    {
        return isFullDel(gainDel) || isPartialDel(gainDel);
    }

    private static boolean isFullDel(@NotNull PurpleGainDeletion gainDel)
    {
        return gainDel.interpretation() == CopyNumberInterpretation.FULL_DEL;
    }

    private static boolean isPartialDel(final @NotNull PurpleGainDeletion gainDel)
    {
        return gainDel.interpretation() == CopyNumberInterpretation.PARTIAL_DEL;
    }

    @NotNull
    private static List<PurpleGeneCopyNumber> mergeSuspectGeneCopyNumberWithLOH(@NotNull PurpleRecord purple)
    {
        List<PurpleLossOfHeterozygosity> reportableGermlineLossOfHeterozygosities = purple.reportableGermlineLossOfHeterozygosities();
        if(reportableGermlineLossOfHeterozygosities == null)
        {
            return purple.suspectGeneCopyNumbersWithLOH();
        }

        Set<String> genesWithReportableGermlineLOH =
                reportableGermlineLossOfHeterozygosities.stream().map(PurpleLossOfHeterozygosity::gene).collect(Collectors.toSet());

        List<PurpleGeneCopyNumber> mergedSuspectGeneCopyNumberWithLOH = Lists.newArrayList();
        for(PurpleGeneCopyNumber somaticSuspectLOH : purple.suspectGeneCopyNumbersWithLOH())
        {
            if(genesWithReportableGermlineLOH.contains(somaticSuspectLOH.gene()))
            {
                PurpleGeneCopyNumber adjustedSuspectLOH =
                        correctForGermlineImpact(somaticSuspectLOH, reportableGermlineLossOfHeterozygosities);
                mergedSuspectGeneCopyNumberWithLOH.add(adjustedSuspectLOH);
            }
            else
            {
                mergedSuspectGeneCopyNumberWithLOH.add(somaticSuspectLOH);
            }
        }

        for(String gene : genesWithReportableGermlineLOH)
        {
            PurpleGeneCopyNumber candidate = findMatchingGeneCopyNumber(gene, purple.allSomaticGeneCopyNumbers());

            if(candidate == null)
            {
                throw new IllegalStateException("No somatic gene copy number found for gene " + gene);
            }

            if(shouldAddSuspectGeneCopyNumberWithLOH(candidate, purple.suspectGeneCopyNumbersWithLOH(), purple.reportableGermlineFullDels()))
            {
                PurpleGeneCopyNumber adjustedGeneCopyNumber = correctForGermlineImpact(candidate, reportableGermlineLossOfHeterozygosities);
                mergedSuspectGeneCopyNumberWithLOH.add(adjustedGeneCopyNumber);
            }
        }
        return mergedSuspectGeneCopyNumberWithLOH;
    }

    private static boolean shouldAddSuspectGeneCopyNumberWithLOH(@NotNull PurpleGeneCopyNumber candidateSuspectGeneCopyNumber,
            @NotNull List<PurpleGeneCopyNumber> suspectGeneCopyNumbersWithLOH, @Nullable List<PurpleGainDeletion> reportableGermlineFullDels)
    {
        String gene = candidateSuspectGeneCopyNumber.gene();

        boolean alreadySuspectGeneCopyNumberWithLOH = findMatchingGeneCopyNumber(gene, suspectGeneCopyNumbersWithLOH) != null;
        boolean hasReportableGermlineFullDels =
                reportableGermlineFullDels != null && reportableGermlineFullDels.stream().anyMatch(l -> l.gene().equals(gene));
        boolean fullyLostInTumor = candidateSuspectGeneCopyNumber.minCopyNumber() < 0.5;

        return !alreadySuspectGeneCopyNumberWithLOH && !hasReportableGermlineFullDels && !fullyLostInTumor;
    }

    @Nullable
    private static PurpleGeneCopyNumber findMatchingGeneCopyNumber(@NotNull String gene,
            @NotNull List<PurpleGeneCopyNumber> geneCopyNumbers)
    {
        for(PurpleGeneCopyNumber geneCopyNumber : geneCopyNumbers)
        {
            if(geneCopyNumber.gene().equals(gene))
            {
                return geneCopyNumber;
            }
        }
        return null;
    }

    @NotNull
    private static PurpleGeneCopyNumber correctForGermlineImpact(@NotNull PurpleGeneCopyNumber geneCopyNumber,
            @NotNull List<PurpleLossOfHeterozygosity> reportableGermlineLossOfHeterozygosities)
    {
        double adjustedMinCopyNumber = adjustMinCopyNumberForGermlineImpact(geneCopyNumber, reportableGermlineLossOfHeterozygosities);
        double adjustedMaxCopyNumber = adjustMaxCopyNumberForGermlineImpact(geneCopyNumber, reportableGermlineLossOfHeterozygosities);
        return ImmutablePurpleGeneCopyNumber.builder()
                .from(geneCopyNumber)
                .minCopyNumber(adjustedMinCopyNumber)
                .minMinorAlleleCopyNumber(0D)
                .maxCopyNumber(adjustedMaxCopyNumber)
                .build();
    }

    private static double adjustMinCopyNumberForGermlineImpact(@NotNull PurpleGeneCopyNumber geneCopyNumber,
            @NotNull List<PurpleLossOfHeterozygosity> reportableGermlineLossOfHeterozygosities)
    {
        OptionalDouble minimumTumorCopyNumber = reportableGermlineLossOfHeterozygosities.stream()
                .filter(d -> d.gene().equals(geneCopyNumber.gene()))
                .mapToDouble(PurpleLossOfHeterozygosity::minCopies)
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

    private static double adjustMaxCopyNumberForGermlineImpact(@NotNull PurpleGeneCopyNumber geneCopyNumber,
            @NotNull List<PurpleLossOfHeterozygosity> reportableGermlineLossOfHeterozygosities)
    {
        OptionalDouble maximumTumorCopyNumber = reportableGermlineLossOfHeterozygosities.stream()
                .filter(d -> d.gene().equals(geneCopyNumber.gene()))
                .mapToDouble(PurpleLossOfHeterozygosity::maxCopies)
                .max();
        if(maximumTumorCopyNumber.isPresent())
        {
            return Math.max(maximumTumorCopyNumber.getAsDouble(), geneCopyNumber.maxCopyNumber());
        }
        else
        {
            return geneCopyNumber.maxCopyNumber();
        }
    }

    @NotNull
    private static PurpleDriver mergeSomaticMutationDriverWithGermline(@NotNull PurpleDriver somaticDriver,
            @NotNull PurpleDriver germlineDriver)
    {
        return ImmutablePurpleDriver.builder()
                .from(somaticDriver)
                .likelihoodMethod(PurpleLikelihoodMethod.HOTSPOT)
                .driverLikelihood(Math.max(somaticDriver.driverLikelihood(), germlineDriver.driverLikelihood()))
                .build();
    }

    @Nullable
    private static PurpleDriver findMatchingGermlineDriver(@NotNull PurpleDriver somaticDriver,
            @Nullable List<PurpleDriver> germlineDrivers)
    {
        if(germlineDrivers == null)
        {
            return null;
        }

        return find(germlineDrivers, PurpleDriverType.GERMLINE_MUTATION, somaticDriver.gene(), somaticDriver.transcript());
    }

    @Nullable
    private static PurpleDriver findMatchingSomaticDriver(@NotNull PurpleDriver germlineDriver, @NotNull List<PurpleDriver> somaticDrivers,
            @NotNull PurpleDriverType somaticTypeToFind)
    {
        return find(somaticDrivers, somaticTypeToFind, germlineDriver.gene(), germlineDriver.transcript());
    }

    @Nullable
    private static PurpleDriver find(@NotNull List<PurpleDriver> drivers, @NotNull PurpleDriverType PurpleDriverTypeToFind,
            @NotNull String geneToFind, @NotNull String transcriptToFind)
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
    private static PurpleDriver convertToSomaticMutationDriver(@NotNull PurpleDriver germlineDriver)
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
    private static PurpleDriver convertToSomaticDeletionDriver(@NotNull PurpleDriver germlineDriver)
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
    static LinxRecord convertLinxGermline(boolean containsTumorCells, @NotNull LinxRecord linx)
    {
        List<LinxSvAnnotation> additionalStructuralVariants = Lists.newArrayList();
        List<LinxBreakend> additionalReportableBreakends = Lists.newArrayList();
        List<LinxHomozygousDisruption> additionalHomozygousDisruptions = Lists.newArrayList();

        if(containsTumorCells)
        {
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
            @NotNull Map<Integer, Integer> svIdMapping, @NotNull Map<Integer, Integer> clusterIdMapping)
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
            @NotNull Map<Integer, Integer> breakendIdMapping, @NotNull Map<Integer, Integer> svIdMapping)
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
    private static Map<Integer, Integer> buildSvIdMapping(@NotNull List<LinxSvAnnotation> allSomaticStructuralVariants,
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
    static int findMaxSvId(@NotNull List<LinxSvAnnotation> structuralVariants)
    {
        int maxSvId = 0;
        for(LinxSvAnnotation structuralVariant : structuralVariants)
        {
            maxSvId = Math.max(maxSvId, structuralVariant.svId());
        }
        return maxSvId;
    }

    @NotNull
    private static Map<Integer, Integer> buildClusterIdMapping(@NotNull List<LinxSvAnnotation> allSomaticStructuralVariants,
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
    static int findMaxClusterId(@NotNull List<LinxSvAnnotation> structuralVariants)
    {
        int maxClusterId = 0;
        for(LinxSvAnnotation structuralVariant : structuralVariants)
        {
            maxClusterId = Math.max(maxClusterId, structuralVariant.clusterId());
        }
        return maxClusterId;
    }

    @NotNull
    private static Map<Integer, Integer> buildBreakendIdMapping(@NotNull List<LinxBreakend> allSomaticBreakends,
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
    static int findMaxBreakendId(@NotNull List<LinxBreakend> breakends)
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
