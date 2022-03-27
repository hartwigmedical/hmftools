package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.groupByImpact;
import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.isReportable;

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
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

import org.apache.logging.log4j.util.Strings;

public class TsgDrivers
{
    private final ReportablePredicate mReportablePredicate;
    private final Map<String, DndsDriverGeneLikelihood> mLikelihoodsByGene;
    private final List<SomaticVariant> mReportableVariants;

    public TsgDrivers(final DriverGenePanel genePanel)
    {
        mLikelihoodsByGene = genePanel.tsgLikelihood();
        mReportablePredicate = new ReportablePredicate(DriverCategory.TSG, genePanel.driverGenes());
        mReportableVariants = Lists.newArrayList();
    }

    public boolean checkVariant(final SomaticVariant variant)
    {
        if(isReportable(mReportablePredicate, variant))
        {
            mReportableVariants.add(variant);
            return true;
        }

        return false;
    }

    public List<DriverCatalog> findDrivers(
            final Map<String,List<GeneCopyNumber>> geneCopyNumberMap,
            final Map<VariantType,Integer> variantTypeCounts, final Map<VariantType,Integer> variantTypeCountsBiallelic)
    {
        final List<DriverCatalog> driverCatalog = Lists.newArrayList();

        final Map<String,List<SomaticVariant>> codingVariants = mReportableVariants.stream()
                .collect(Collectors.groupingBy(SomaticVariant::gene));

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

        SomaticVariant topVariant = geneVariants.get(0);

        Map<DriverImpact,Integer> variantCounts = groupByImpact(geneVariants);
        int missenseVariants = variantCounts.getOrDefault(DriverImpact.MISSENSE, 0);
        int nonsenseVariants = variantCounts.getOrDefault(DriverImpact.NONSENSE, 0);
        int spliceVariants = variantCounts.getOrDefault(DriverImpact.SPLICE, 0);
        int inframeVariants = variantCounts.getOrDefault(DriverImpact.INFRAME, 0);
        int frameshiftVariants = variantCounts.getOrDefault(DriverImpact.FRAMESHIFT, 0);

        final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .chromosome(topVariant.chromosome())
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

        if(geneVariants.stream().anyMatch(x -> x.biallelic() && !DriverImpact.isMissense(x.type(), x.variantImpact().CanonicalCodingEffect)))
        {
            return builder.likelihoodMethod(LikelihoodMethod.BIALLELIC).build();
        }

        final DriverImpact firstImpact = DriverImpact.select(topVariant.type(), topVariant.variantImpact().CanonicalCodingEffect);
        final DndsDriverImpactLikelihood firstImpactLikelihood = likelihood.select(firstImpact);
        final int firstVariantTypeCount = variantCount(topVariant.biallelic(), topVariant, standardCounts, biallelicCounts);

        int nonBiallelicMissenseCount = standardCounts.getOrDefault(VariantType.SNP, 0);

        if(geneVariants.size() == 1)
        {
            double singleHit = singleHit(firstVariantTypeCount, firstImpactLikelihood);
            double substituteFirst =
                    firstImpact != DriverImpact.MISSENSE ? singleHit(nonBiallelicMissenseCount, likelihood.missense()) : singleHit;

            return builder.driverLikelihood(Math.max(singleHit, substituteFirst)).build();
        }

        // MultiHit
        SomaticVariant secondVariant = geneVariants.get(1);

        final DriverImpact secondImpact = DriverImpact.select(secondVariant.type(), secondVariant.variantImpact().CanonicalCodingEffect);

        final DndsDriverImpactLikelihood secondImpactLikelihood = likelihood.select(secondImpact);

        final int secondVariantTypeCount =
                variantCount(secondVariant.biallelic(), secondVariant, standardCounts, biallelicCounts);

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

    private static int variantCount(
            boolean useBiallelic, final SomaticVariant variant,
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

        return variant.type() == VariantType.INDEL ?
                map.getOrDefault(VariantType.INDEL, 0) : map.getOrDefault(VariantType.SNP, 0);
    }
}
