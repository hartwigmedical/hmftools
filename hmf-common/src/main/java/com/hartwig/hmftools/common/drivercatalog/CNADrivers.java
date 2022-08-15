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
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;

public class CNADrivers
{
    private static final double MIN_COPY_NUMBER_RELATIVE_INCREASE = 3;

    public static final double MAX_COPY_NUMBER_DEL = 0.5;

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

    public static int deletedGenes(final List<GeneCopyNumber> geneCopyNumbers)
    {
        return (int) geneCopyNumbers.stream()
                .filter(x -> !HumanChromosome.fromString(x.chromosome()).equals(HumanChromosome._Y)
                        && Doubles.lessThan(x.minCopyNumber(), MAX_COPY_NUMBER_DEL))
                .count();
    }

    public List<DriverCatalog> amplifications(final double ploidy, final List<GeneCopyNumber> geneCopyNumbers)
    {
        return amplifications(ploidy, geneCopyNumbers, true);
    }

    public List<DriverCatalog> amplifications(final double ploidy, final List<GeneCopyNumber> geneCopyNumbers, boolean allowPartials)
    {
        Predicate<GeneCopyNumber> targetPredicate = x -> mAmplificationTargets.containsKey(x.geneName());
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
                    result.add(createAmpDriver(geneCopyNumber));
                }
                else if(allowPartials && maxCopyNumberPredicate.test(geneCopyNumber))
                {
                    result.add(createPartialAmpDriver(geneCopyNumber));
                }
            }
        }

        return result;
    }

    public List<DriverCatalog> deletions(final List<GeneCopyNumber> geneCopyNumbers)
    {
        List<DriverCatalog> drivers = Lists.newArrayList();

        boolean checkQcStatus = mQcStatus.contains(PurpleQCStatus.WARN_DELETED_GENES)
                || mQcStatus.contains(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE);

        for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
        {
            if(geneCopyNumber.minCopyNumber() >= MAX_COPY_NUMBER_DEL)
                continue;

            if(!mDeletionTargets.containsKey(geneCopyNumber.geneName()))
                continue;

            if(checkQcStatus)
            {
                if(!supportedByTwoSVs(geneCopyNumber) && !shortAndSupportedByOneSVAndMere(geneCopyNumber))
                    continue;
            }

            drivers.add(createDelDriver(geneCopyNumber));
        }

        return drivers;
    }

    static boolean supportedByOneSV(final GeneCopyNumber geneCopyNumber)
    {
        return geneCopyNumber.minRegionStartSupport().isSV() || geneCopyNumber.minRegionEndSupport().isSV();
    }

    static boolean supportedByTwoSVs(final GeneCopyNumber geneCopyNumber)
    {
        return geneCopyNumber.minRegionStartSupport().isSV() && geneCopyNumber.minRegionEndSupport().isSV();
    }

    static boolean shortAndSupportedByOneSVAndMere(final GeneCopyNumber geneCopyNumber)
    {
        if(geneCopyNumber.minRegionBases() >= 10_000_000)
            return false;

        if(MERE.contains(geneCopyNumber.minRegionStartSupport()))
            return geneCopyNumber.minRegionEndSupport().isSV();

        if(MERE.contains(geneCopyNumber.minRegionEndSupport()))
            return geneCopyNumber.minRegionStartSupport().isSV();

        return false;
    }

    private DriverCatalog createAmpDriver(final GeneCopyNumber geneCopyNumber)
    {
        DriverGene driverGene = mAmplificationTargets.get(geneCopyNumber.geneName());
        return cnaDriver(driverGene.likelihoodType(), DriverType.AMP, LikelihoodMethod.AMP, false, geneCopyNumber);
    }

    private DriverCatalog createPartialAmpDriver(final GeneCopyNumber geneCopyNumber)
    {
        DriverGene driverGene = mAmplificationTargets.get(geneCopyNumber.geneName());
        return cnaDriver(driverGene.likelihoodType(), DriverType.PARTIAL_AMP, LikelihoodMethod.AMP, false, geneCopyNumber);
    }

    private DriverCatalog createDelDriver(final GeneCopyNumber geneCopyNumber)
    {
        DriverGene driverGene = mDeletionTargets.get(geneCopyNumber.geneName());
        return cnaDriver(driverGene.likelihoodType(), DriverType.DEL, LikelihoodMethod.DEL, true, geneCopyNumber);
    }

    private DriverCatalog cnaDriver(
            DriverCategory category, DriverType driver, LikelihoodMethod likelihoodMethod, boolean biallelic,
            final GeneCopyNumber geneCopyNumber)
    {
        return ImmutableDriverCatalog.builder()
                .chromosome(geneCopyNumber.chromosome())
                .chromosomeBand(geneCopyNumber.chromosomeBand())
                .gene(geneCopyNumber.geneName())
                .transcript(geneCopyNumber.transName())
                .isCanonical(geneCopyNumber.isCanonical())
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
                .minCopyNumber(geneCopyNumber.minCopyNumber())
                .maxCopyNumber(geneCopyNumber.maxCopyNumber())
                .build();
    }
}
