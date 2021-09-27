package com.hartwig.hmftools.purple.drivers;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFactory;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverImpact;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.region.TranscriptRegion;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class TsgDrivers
{
    private final ReportablePredicate mReportablePredicate;
    private final Map<String, DndsDriverGeneLikelihood> mLikelihoodsByGene;

    public TsgDrivers(@NotNull DriverGenePanel genePanel)
    {
        mLikelihoodsByGene = genePanel.tsgLikelihood();
        mReportablePredicate = new ReportablePredicate(DriverCategory.TSG, genePanel);
    }

    @NotNull
    List<DriverCatalog> drivers(@NotNull final List<SomaticVariant> variants, @NotNull final List<GeneCopyNumber> geneCopyNumberList,
            final Map<VariantType, Long> variantTypeCounts, final Map<VariantType, Long> variantTypeCountsBiallelic)
    {
        final Map<String, GeneCopyNumber> geneCopyNumbers =
                geneCopyNumberList.stream().collect(Collectors.toMap(TranscriptRegion::gene, x -> x));
        final List<DriverCatalog> driverCatalog = Lists.newArrayList();

        final Map<String, List<SomaticVariant>> codingVariants = codingVariantsByGene(variants);

        for(String gene : codingVariants.keySet())
        {
            final DndsDriverGeneLikelihood likelihood = mLikelihoodsByGene.get(gene);

            final List<SomaticVariant> geneVariants = codingVariants.get(gene);
            driverCatalog.add(geneDriver(likelihood,
                    geneVariants,
                    variantTypeCounts,
                    variantTypeCountsBiallelic,
                    geneCopyNumbers.get(gene)));
        }

        return driverCatalog;
    }

    @NotNull
    static DriverCatalog geneDriver(@NotNull final DndsDriverGeneLikelihood likelihood, @NotNull final List<SomaticVariant> geneVariants,
            @NotNull final Map<VariantType, Long> standardCounts, @NotNull final Map<VariantType, Long> biallelicCounts,
            @Nullable GeneCopyNumber geneCopyNumber)
    {
        geneVariants.sort(new TsgImpactComparator());

        final Map<DriverImpact, Long> variantCounts = DriverCatalogFactory.driverImpactCount(geneVariants);
        long missenseVariants = variantCounts.getOrDefault(DriverImpact.MISSENSE, 0L);
        long nonsenseVariants = variantCounts.getOrDefault(DriverImpact.NONSENSE, 0L);
        long spliceVariants = variantCounts.getOrDefault(DriverImpact.SPLICE, 0L);
        long inframeVariants = variantCounts.getOrDefault(DriverImpact.INFRAME, 0L);
        long frameshiftVariants = variantCounts.getOrDefault(DriverImpact.FRAMESHIFT, 0L);

        final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .chromosome(geneVariants.get(0).chromosome())
                .chromosomeBand(geneCopyNumber == null ? Strings.EMPTY : geneCopyNumber.chromosomeBand())
                .gene(likelihood.gene())
                .driver(DriverType.MUTATION)
                .category(DriverCategory.TSG)
                .driverLikelihood(1)
                .missense(missenseVariants)
                .nonsense(nonsenseVariants)
                .splice(spliceVariants)
                .inframe(inframeVariants)
                .frameshift(frameshiftVariants)
                .biallelic(geneVariants.stream().anyMatch(SomaticVariant::biallelic))
                .minCopyNumber(geneCopyNumber == null ? 0 : geneCopyNumber.minCopyNumber())
                .maxCopyNumber(geneCopyNumber == null ? 0 : geneCopyNumber.maxCopyNumber())
                .likelihoodMethod(LikelihoodMethod.DNDS);

        if(geneVariants.stream().anyMatch(SomaticVariant::isHotspot))
        {
            return builder.likelihoodMethod(LikelihoodMethod.HOTSPOT).build();
        }

        if(geneVariants.stream().anyMatch(x -> x.biallelic() && !DriverImpact.isMissense(x)))
        {
            return builder.likelihoodMethod(LikelihoodMethod.BIALLELIC).build();
        }

        final DriverImpact firstImpact = DriverImpact.select(geneVariants.get(0));
        final DndsDriverImpactLikelihood firstImpactLikelihood = likelihood.select(firstImpact);
        final long firstVariantTypeCount =
                variantCount(geneVariants.get(0).biallelic(), geneVariants.get(0), standardCounts, biallelicCounts);

        long nonBiallelicMissenseCount = standardCounts.getOrDefault(VariantType.SNP, 0L);

        if(geneVariants.size() == 1)
        {
            double singleHit = singleHit(firstVariantTypeCount, firstImpactLikelihood);
            double substituteFirst =
                    firstImpact != DriverImpact.MISSENSE ? singleHit(nonBiallelicMissenseCount, likelihood.missense()) : singleHit;

            return builder.driverLikelihood(Math.max(singleHit, substituteFirst)).build();
        }

        // MultiHit
        final DriverImpact secondImpact = DriverImpact.select(geneVariants.get(1));
        final DndsDriverImpactLikelihood secondImpactLikelihood = likelihood.select(secondImpact);
        final long secondVariantTypeCount =
                variantCount(geneVariants.get(1).biallelic(), geneVariants.get(1), standardCounts, biallelicCounts);
        double multiHit = multiHit(firstVariantTypeCount, secondVariantTypeCount, firstImpactLikelihood, secondImpactLikelihood);

        double substituteFirst = firstImpact == DriverImpact.MISSENSE
                ? multiHit
                : multiHit(nonBiallelicMissenseCount, secondVariantTypeCount, likelihood.missense(), secondImpactLikelihood);

        double substituteSecond = secondImpact == DriverImpact.MISSENSE
                ? multiHit
                : multiHit(firstVariantTypeCount, nonBiallelicMissenseCount, firstImpactLikelihood, likelihood.missense());

        double substituteBoth = firstImpact == DriverImpact.MISSENSE || secondImpact == DriverImpact.MISSENSE
                ? multiHit
                : multiHit(nonBiallelicMissenseCount, nonBiallelicMissenseCount, likelihood.missense(), likelihood.missense());

        double combinedResult = Math.max(Math.max(substituteFirst, substituteSecond), substituteBoth);
        return builder.driverLikelihood(combinedResult).build();
    }

    private static double multiHit(long firstVariantTypeCount, long secondVariantTypeCount,
            @NotNull final DndsDriverImpactLikelihood firstLikelihood, @NotNull final DndsDriverImpactLikelihood secondLikelihood)
    {
        return DriverCatalogFactory.probabilityDriverVariant(firstVariantTypeCount,
                secondVariantTypeCount,
                firstLikelihood,
                secondLikelihood);
    }

    private static double singleHit(long sampleCount, @NotNull final DndsDriverImpactLikelihood likelihood)
    {
        return DriverCatalogFactory.probabilityDriverVariant(sampleCount, likelihood);
    }

    @NotNull
    private <T extends SomaticVariant> Map<String, List<T>> codingVariantsByGene(@NotNull final List<T> variants)
    {
        return variants.stream().filter(mReportablePredicate).collect(Collectors.groupingBy(SomaticVariant::gene));
    }

    private static long variantCount(boolean useBiallelic, @NotNull final SomaticVariant variant,
            @NotNull final Map<VariantType, Long> standard, @NotNull final Map<VariantType, Long> biallelic)
    {
        final Map<VariantType, Long> map;
        if(!useBiallelic)
        {
            map = standard;
        }
        else if(variant.biallelic())
        {
            map = biallelic;
        }
        else
        {
            map = standard;
        }

        return variant.type() == VariantType.INDEL ? map.getOrDefault(VariantType.INDEL, 0L) : map.getOrDefault(VariantType.SNP, 0L);
    }
}
