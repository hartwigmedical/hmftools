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
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleDriver;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleFit;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleFittedPurityMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleLikelihoodMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleLossOfHeterozygosity;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;

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
        List<PurpleDriver> mergedDrivers;
        List<PurpleVariant> additionalReportableVariants;
        List<PurpleGainLoss> mergedAllSomaticGainsLosses;
        List<PurpleGainLoss> mergedReportableSomaticGainsLosses;
        List<PurpleGeneCopyNumber> mergedSuspectGeneCopyNumbersWithLOH;
        if(containsTumorCells)
        {
            mergedDrivers =
                    mergeGermlineDriversIntoSomatic(purple.somaticDrivers(), purple.germlineDrivers(), purple.reportableGermlineFullLosses());
            additionalReportableVariants = toSomaticVariants(purple.reportableGermlineVariants());
            mergedAllSomaticGainsLosses = mergeGermlineFullLosses(purple.reportableGermlineFullLosses(), purple.allSomaticGainsLosses());
            mergedReportableSomaticGainsLosses =
                    mergeGermlineFullLosses(purple.reportableGermlineFullLosses(), purple.reportableSomaticGainsLosses());
            mergedSuspectGeneCopyNumbersWithLOH = mergeSuspectGeneCopyNumberWithLOH(purple);
        }
        else
        {
            mergedDrivers = purple.somaticDrivers();
            additionalReportableVariants = Lists.newArrayList();
            mergedAllSomaticGainsLosses = purple.allSomaticGainsLosses();
            mergedReportableSomaticGainsLosses = purple.reportableSomaticGainsLosses();
            mergedSuspectGeneCopyNumbersWithLOH = purple.suspectGeneCopyNumbersWithLOH();
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
                .allSomaticGainsLosses(mergedAllSomaticGainsLosses)
                .reportableSomaticGainsLosses(mergedReportableSomaticGainsLosses)
                .allGermlineDeletions(null)
                .allGermlineFullLosses(null)
                .reportableGermlineFullLosses(null)
                .suspectGeneCopyNumbersWithLOH(mergedSuspectGeneCopyNumbersWithLOH)
                .allGermlineLossOfHeterozygosities(null)
                .reportableGermlineLossOfHeterozygosities(null)
                .build();
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
            @Nullable List<PurpleDriver> germlineDrivers, @Nullable List<PurpleGainLoss> reportableGermlineFullLosses)
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
                        && reportableGermlineFullLosses != null && reportableGermlineFullLosses.stream()
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
    private static List<PurpleGainLoss> mergeGermlineFullLosses(@Nullable List<PurpleGainLoss> reportableGermlineFullLosses,
            @NotNull List<PurpleGainLoss> somaticGainLosses)
    {
        if(reportableGermlineFullLosses == null)
        {
            return somaticGainLosses;
        }

        List<PurpleGainLoss> mergedGainLosses = Lists.newArrayList();
        for(PurpleGainLoss somaticGainLoss : somaticGainLosses)
        {
            mergedGainLosses.add(getGermlineCorrectedGainLoss(somaticGainLoss, reportableGermlineFullLosses));
        }

        for(PurpleGainLoss germlineLoss : reportableGermlineFullLosses)
        {
            boolean hasMatchingSomaticLoss = somaticGainLosses.stream().anyMatch(s -> lossesMatch(germlineLoss, s));
            if(!hasMatchingSomaticLoss)
            {
                mergedGainLosses.add(germlineLoss);
            }
        }

        return mergedGainLosses;
    }

    @NotNull
    private static PurpleGainLoss getGermlineCorrectedGainLoss(@NotNull PurpleGainLoss somaticGainLoss,
            @NotNull List<PurpleGainLoss> reportableGermlineFullLosses)
    {
        if(!isLoss(somaticGainLoss))
        {
            return somaticGainLoss;
        }

        List<PurpleGainLoss> matchingGermlineLosses =
                reportableGermlineFullLosses.stream().filter(l -> lossesMatch(l, somaticGainLoss)).collect(Collectors.toList());
        if(matchingGermlineLosses.isEmpty())
        {
            return somaticGainLoss;
        }
        else if(matchingGermlineLosses.size() == 1)
        {
            return mergeLosses(somaticGainLoss, matchingGermlineLosses.get(0));
        }
        else
        {
            throw new IllegalStateException(
                    "More than one reportable germline loss for transcript " + somaticGainLoss.transcript() + " of gene "
                            + somaticGainLoss.gene());
        }
    }

    private static boolean lossesMatch(@NotNull PurpleGainLoss germlineLoss, @NotNull PurpleGainLoss somaticGainLoss)
    {
        return isLoss(somaticGainLoss) && germlineLoss.gene().equals(somaticGainLoss.gene()) && germlineLoss.transcript()
                .equals(somaticGainLoss.transcript());
    }

    @NotNull
    private static PurpleGainLoss mergeLosses(@NotNull PurpleGainLoss somaticLoss, @NotNull PurpleGainLoss germlineLoss)
    {
        double minCopies = Math.min(somaticLoss.minCopies(), germlineLoss.minCopies());

        CopyNumberInterpretation interpretation;
        double maxCopies;
        if(somaticLoss.interpretation() == germlineLoss.interpretation())
        {
            interpretation = somaticLoss.interpretation();
            maxCopies = Math.max(somaticLoss.maxCopies(), germlineLoss.maxCopies());
        }
        else if(isFullLoss(somaticLoss))
        {
            interpretation = CopyNumberInterpretation.FULL_LOSS;
            maxCopies = somaticLoss.maxCopies();
        }
        else
        {
            // germline full loss and somatic partial loss
            interpretation = CopyNumberInterpretation.FULL_LOSS;
            maxCopies = germlineLoss.maxCopies();
        }

        return ImmutablePurpleGainLoss.builder()
                .interpretation(interpretation)
                .chromosome(somaticLoss.chromosome())
                .chromosomeBand(somaticLoss.chromosomeBand())
                .gene(somaticLoss.gene())
                .transcript(somaticLoss.transcript())
                .isCanonical(somaticLoss.isCanonical())
                .minCopies(minCopies)
                .maxCopies(maxCopies)
                .build();
    }

    private static boolean isLoss(@NotNull PurpleGainLoss gainLoss)
    {
        return isFullLoss(gainLoss) || isPartialLoss(gainLoss);
    }

    private static boolean isFullLoss(@NotNull PurpleGainLoss gainLoss)
    {
        return gainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS;
    }

    private static boolean isPartialLoss(final @NotNull PurpleGainLoss gainLoss)
    {
        return gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS;
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

            if(shouldAddSuspectGeneCopyNumberWithLOH(candidate, purple.suspectGeneCopyNumbersWithLOH(), purple.reportableGermlineFullLosses()))
            {
                PurpleGeneCopyNumber adjustedGeneCopyNumber = correctForGermlineImpact(candidate, reportableGermlineLossOfHeterozygosities);
                mergedSuspectGeneCopyNumberWithLOH.add(adjustedGeneCopyNumber);
            }
        }
        return mergedSuspectGeneCopyNumberWithLOH;
    }

    private static boolean shouldAddSuspectGeneCopyNumberWithLOH(@NotNull PurpleGeneCopyNumber candidateSuspectGeneCopyNumber,
            @NotNull List<PurpleGeneCopyNumber> suspectGeneCopyNumbersWithLOH, @Nullable List<PurpleGainLoss> reportableGermlineFullLosses)
    {
        String gene = candidateSuspectGeneCopyNumber.gene();

        boolean alreadySuspectGeneCopyNumberWithLOH = findMatchingGeneCopyNumber(gene, suspectGeneCopyNumbersWithLOH) != null;
        boolean hasReportableGermlineFullLoss =
                reportableGermlineFullLosses != null && reportableGermlineFullLosses.stream().anyMatch(l -> l.gene().equals(gene));
        boolean fullyLostInTumor = candidateSuspectGeneCopyNumber.minCopyNumber() < 0.5;

        return !alreadySuspectGeneCopyNumberWithLOH && !hasReportableGermlineFullLoss && !fullyLostInTumor;
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
