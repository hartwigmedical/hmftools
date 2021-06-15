package com.hartwig.hmftools.common.drivercatalog;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public class CNADrivers
{
    private static final double MIN_COPY_NUMBER_RELATIVE_INCREASE = 3;
    private static final double MAX_COPY_NUMBER_DEL = 0.5;

    private static final Set<SegmentSupport> MERE = Sets.newHashSet(SegmentSupport.CENTROMERE, SegmentSupport.TELOMERE);

    private final Set<PurpleQCStatus> mQcStatus;
    private final Map<String, DriverGene> mAmplificationTargets;
    private final Map<String, DriverGene> mDeletionTargets;

    public CNADrivers(Set<PurpleQCStatus> qcStatus, DriverGenePanel panel)
    {
        mQcStatus = qcStatus;
        mDeletionTargets = panel.deletionTargets().stream().collect(Collectors.toMap(DriverGene::gene, x -> x));
        mAmplificationTargets = panel.amplificationTargets().stream().collect(Collectors.toMap(DriverGene::gene, x -> x));
    }

    public static int deletedGenes(@NotNull final List<GeneCopyNumber> geneCopyNumbers)
    {
        return (int) geneCopyNumbers.stream()
                .filter(x -> !HumanChromosome.fromString(x.chromosome()).equals(HumanChromosome._Y) && x.germlineHet2HomRegions() == 0
                        && x.germlineHomRegions() == 0 && Doubles.lessThan(x.minCopyNumber(), MAX_COPY_NUMBER_DEL))
                .count();
    }

    @NotNull
    public List<DriverCatalog> amplifications(final double ploidy, @NotNull final List<GeneCopyNumber> geneCopyNumbers)
    {
        Predicate<GeneCopyNumber> targetPredicate = x -> mAmplificationTargets.containsKey(x.gene());
        Predicate<GeneCopyNumber> qcStatusPredicate =
                mQcStatus.contains(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE) ? CNADrivers::supportedByOneSV : x -> true;

        Predicate<GeneCopyNumber> minCopyNumberPredicate = x -> x.minCopyNumber() / ploidy > MIN_COPY_NUMBER_RELATIVE_INCREASE;
        Predicate<GeneCopyNumber> maxCopyNumberPredicate = x -> x.maxCopyNumber() / ploidy > MIN_COPY_NUMBER_RELATIVE_INCREASE;

        List<DriverCatalog> result = Lists.newArrayList();
        for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
        {
            if(qcStatusPredicate.test(geneCopyNumber) && targetPredicate.test(geneCopyNumber))
            {
                if(minCopyNumberPredicate.test(geneCopyNumber))
                {
                    result.add(amp(geneCopyNumber));
                }
                else if(maxCopyNumberPredicate.test(geneCopyNumber))
                {
                    result.add(partialAmp(geneCopyNumber));
                }
            }
        }

        return result;
    }

    @NotNull
    public List<DriverCatalog> deletions(@NotNull final List<GeneCopyNumber> geneCopyNumbers)
    {
        Predicate<GeneCopyNumber> qcStatusPredicate =
                mQcStatus.contains(PurpleQCStatus.WARN_DELETED_GENES) || mQcStatus.contains(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE) ? x ->
                        supportedByTwoSVs(x) || shortAndSupportedByOneSVAndMere(x) : x -> true;

        return geneCopyNumbers.stream()
                .filter(x -> x.minCopyNumber() < MAX_COPY_NUMBER_DEL)
                .filter(x -> mDeletionTargets.containsKey(x.gene()))
                .filter(x -> x.germlineHet2HomRegions() == 0 && x.germlineHomRegions() == 0)
                .filter(qcStatusPredicate)
                .map(this::del)
                .collect(Collectors.toList());
    }

    static boolean supportedByOneSV(@NotNull final GeneCopyNumber geneCopyNumber)
    {
        return geneCopyNumber.minRegionStartSupport().isSV() || geneCopyNumber.minRegionEndSupport().isSV();
    }

    static boolean supportedByTwoSVs(@NotNull final GeneCopyNumber geneCopyNumber)
    {
        return geneCopyNumber.minRegionStartSupport().isSV() && geneCopyNumber.minRegionEndSupport().isSV();
    }

    static boolean shortAndSupportedByOneSVAndMere(@NotNull final GeneCopyNumber geneCopyNumber)
    {
        if(geneCopyNumber.minRegionBases() >= 10_000_000)
        {
            return false;
        }

        if(MERE.contains(geneCopyNumber.minRegionStartSupport()))
        {
            return geneCopyNumber.minRegionEndSupport().isSV();
        }

        if(MERE.contains(geneCopyNumber.minRegionEndSupport()))
        {
            return geneCopyNumber.minRegionStartSupport().isSV();
        }

        return false;
    }

    @NotNull
    private DriverCatalog amp(GeneCopyNumber x)
    {
        DriverGene driverGene = mAmplificationTargets.get(x.gene());
        return cnaDriver(driverGene.likelihoodType(), DriverType.AMP, LikelihoodMethod.AMP, false, x);
    }

    @NotNull
    private DriverCatalog partialAmp(GeneCopyNumber x)
    {
        DriverGene driverGene = mAmplificationTargets.get(x.gene());
        return cnaDriver(driverGene.likelihoodType(), DriverType.PARTIAL_AMP, LikelihoodMethod.AMP, false, x);
    }

    @NotNull
    private DriverCatalog del(GeneCopyNumber x)
    {
        DriverGene driverGene = mDeletionTargets.get(x.gene());
        return cnaDriver(driverGene.likelihoodType(), DriverType.DEL, LikelihoodMethod.DEL, true, x);
    }

    @NotNull
    private DriverCatalog cnaDriver(DriverCategory category, DriverType driver, LikelihoodMethod likelihoodMethod, boolean biallelic,
            GeneCopyNumber x)
    {
        return ImmutableDriverCatalog.builder()
                .chromosome(x.chromosome())
                .chromosomeBand(x.chromosomeBand())
                .gene(x.gene())
                .missense(0)
                .nonsense(0)
                .inframe(0)
                .frameshift(0)
                .splice(0)
                .driverLikelihood(1)
                .driver(driver)
                .likelihoodMethod(likelihoodMethod)
                .category(category)
                .biallelic(biallelic)
                .minCopyNumber(x.minCopyNumber())
                .maxCopyNumber(x.maxCopyNumber())
                .build();
    }
}
