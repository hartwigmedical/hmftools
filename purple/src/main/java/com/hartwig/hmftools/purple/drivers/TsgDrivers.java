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
import com.hartwig.hmftools.common.drivercatalog.panel.ReportablePredicate;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.util.Strings;

public class TsgDrivers
{
    private final ReportablePredicate mReportablePredicate;
    private final Map<String, DndsDriverGeneLikelihood> mLikelihoodsByGene;

    public TsgDrivers(final DriverGenePanel genePanel)
    {
        mLikelihoodsByGene = genePanel.tsgLikelihood();
        mReportablePredicate = new ReportablePredicate(DriverCategory.TSG, genePanel.driverGenes());
    }

    public List<DriverCatalog> drivers(
            final List<SomaticVariant> variants, final Map<String,List<GeneCopyNumber>> geneCopyNumberMap,
            final Map<VariantType,Integer> variantTypeCounts, final Map<VariantType,Integer> variantTypeCountsBiallelic)
    {
        final List<DriverCatalog> driverCatalog = Lists.newArrayList();

        final Map<String, List<SomaticVariant>> codingVariants = codingVariantsByGene(variants);

        for(String gene : codingVariants.keySet())
        {
            final DndsDriverGeneLikelihood likelihood = mLikelihoodsByGene.get(gene);

            final List<SomaticVariant> geneVariants = codingVariants.get(gene);

            List<GeneCopyNumber> geneCopyNumbers = geneCopyNumberMap.get(gene);

            if(geneCopyNumbers == null)
                continue;

            for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
            {
                driverCatalog.add(geneDriver(likelihood, geneVariants, variantTypeCounts, variantTypeCountsBiallelic, geneCopyNumber));
            }
        }

        return driverCatalog;
    }

    public static DriverCatalog geneDriver(
            final DndsDriverGeneLikelihood likelihood, final List<SomaticVariant> geneVariants, final Map<VariantType, Integer> standardCounts,
            final Map<VariantType,Integer> biallelicCounts, final GeneCopyNumber geneCopyNumber)
    {
        geneVariants.sort(new TsgImpactComparator());

        final Map<DriverImpact,Integer> variantCounts = DriverCatalogFactory.driverImpactCount(geneVariants);
        int missenseVariants = variantCounts.getOrDefault(DriverImpact.MISSENSE, 0);
        int nonsenseVariants = variantCounts.getOrDefault(DriverImpact.NONSENSE, 0);
        int spliceVariants = variantCounts.getOrDefault(DriverImpact.SPLICE, 0);
        int inframeVariants = variantCounts.getOrDefault(DriverImpact.INFRAME, 0);
        int frameshiftVariants = variantCounts.getOrDefault(DriverImpact.FRAMESHIFT, 0);

        final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .chromosome(geneVariants.get(0).chromosome())
                .chromosomeBand(geneCopyNumber == null ? Strings.EMPTY : geneCopyNumber.chromosomeBand())
                .gene(likelihood.gene())
                .transcript(geneCopyNumber != null ? geneCopyNumber.transName() : "")
                .isCanonical(geneCopyNumber != null ? geneCopyNumber.isCanonical() : true)
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
        final int firstVariantTypeCount =
                variantCount(geneVariants.get(0).biallelic(), geneVariants.get(0), standardCounts, biallelicCounts);

        int nonBiallelicMissenseCount = standardCounts.getOrDefault(VariantType.SNP, 0);

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
        final int secondVariantTypeCount =
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

    private static double multiHit(int firstVariantTypeCount, int secondVariantTypeCount,
            final DndsDriverImpactLikelihood firstLikelihood, final DndsDriverImpactLikelihood secondLikelihood)
    {
        return DriverCatalogFactory.probabilityDriverVariant(firstVariantTypeCount,
                secondVariantTypeCount,
                firstLikelihood,
                secondLikelihood);
    }

    private static double singleHit(int sampleCount, final DndsDriverImpactLikelihood likelihood)
    {
        return DriverCatalogFactory.probabilityDriverVariant(sampleCount, likelihood);
    }

    private <T extends SomaticVariant> Map<String, List<T>> codingVariantsByGene(final List<T> variants)
    {
        return variants.stream().filter(mReportablePredicate).collect(Collectors.groupingBy(SomaticVariant::gene));
    }

    private static int variantCount(boolean useBiallelic, final SomaticVariant variant,
            final Map<VariantType, Integer> standard, final Map<VariantType, Integer> biallelic)
    {
        final Map<VariantType, Integer> map;
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

        return variant.type() == VariantType.INDEL ? map.getOrDefault(VariantType.INDEL, 0) : map.getOrDefault(VariantType.SNP, 0);
    }
}
