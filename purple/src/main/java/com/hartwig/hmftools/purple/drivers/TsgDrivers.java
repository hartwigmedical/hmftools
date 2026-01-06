package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.purple.drivers.DndsCalculator.probabilityDriverVariant;
import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.getWorstReportableCodingEffect;
import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.groupByImpact;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.DriverImpact;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.driver.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.driver.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.purple.DriverGeneResource;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

public class TsgDrivers extends SomaticVariantDriverFinder
{
    public TsgDrivers(final DriverGeneResource genePanel)
    {
        super(genePanel, DriverCategory.TSG);
    }

    public DriverCatalog createDriverCatalog(
            final List<SomaticVariant> geneVariants, final Map<VariantType,Integer> variantTypeCounts,
            final Map<VariantType,Integer> biallelicCounts, final GeneCopyNumber geneCopyNumber,
            final DndsDriverGeneLikelihood dndsLikelihood)
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
                .chromosomeBand(geneCopyNumber.ChromosomeBand)
                .gene(geneCopyNumber.geneName())
                .transcript(geneCopyNumber.TransName)
                .isCanonical(geneCopyNumber.IsCanonical)
                .driver(DriverType.MUTATION)
                .category(DriverCategory.TSG)
                .driverLikelihood(1)
                .missense(missenseVariants)
                .nonsense(nonsenseVariants)
                .splice(spliceVariants)
                .inframe(inframeVariants)
                .frameshift(frameshiftVariants)
                .biallelic(geneVariants.stream().anyMatch(SomaticVariant::biallelic))
                .minCopyNumber(geneCopyNumber.minCopyNumber())
                .maxCopyNumber(geneCopyNumber.maxCopyNumber())
                .likelihoodMethod(LikelihoodMethod.DNDS)
                .reportedStatus(ReportedStatus.REPORTED);

        if(geneVariants.stream().anyMatch(SomaticVariant::isHotspot))
        {
            return builder.likelihoodMethod(LikelihoodMethod.HOTSPOT).build();
        }

        if(geneVariants.stream().anyMatch(x -> x.biallelic() && !DriverImpact.isMissense(x.type(), x.variantImpact().CanonicalCodingEffect)))
        {
            return builder.likelihoodMethod(LikelihoodMethod.BIALLELIC).build();
        }

        CodingEffect firstCodingEffect = getWorstReportableCodingEffect(topVariant.variantImpact());
        final DriverImpact firstImpact = DriverImpact.select(topVariant.type(), firstCodingEffect);
        final DndsDriverImpactLikelihood firstImpactLikelihood = dndsLikelihood.select(firstImpact);
        final int firstVariantTypeCount = variantCount(topVariant.biallelic(), topVariant, variantTypeCounts, biallelicCounts);

        int nonBiallelicMissenseCount = variantTypeCounts.getOrDefault(VariantType.SNP, 0);

        if(geneVariants.size() == 1)
        {
            double singleHit = singleHit(firstVariantTypeCount, firstImpactLikelihood);
            double substituteFirst =
                    firstImpact != DriverImpact.MISSENSE ? singleHit(nonBiallelicMissenseCount, dndsLikelihood.missense()) : singleHit;

            return builder.driverLikelihood(Math.max(singleHit, substituteFirst)).build();
        }

        // MultiHit
        SomaticVariant secondVariant = geneVariants.get(1);

        CodingEffect secondCodingEffect = getWorstReportableCodingEffect(secondVariant.variantImpact());
        DriverImpact secondImpact = DriverImpact.select(secondVariant.type(), secondCodingEffect);

        final DndsDriverImpactLikelihood secondImpactLikelihood = dndsLikelihood.select(secondImpact);

        final int secondVariantTypeCount = variantCount(secondVariant.biallelic(), secondVariant, variantTypeCounts, biallelicCounts);

        double multiHit = multiHit(firstVariantTypeCount, secondVariantTypeCount, firstImpactLikelihood, secondImpactLikelihood);

        double substituteFirst = firstImpact == DriverImpact.MISSENSE
                ? multiHit
                : multiHit(nonBiallelicMissenseCount, secondVariantTypeCount, dndsLikelihood.missense(), secondImpactLikelihood);

        double substituteSecond = secondImpact == DriverImpact.MISSENSE
                ? multiHit
                : multiHit(firstVariantTypeCount, nonBiallelicMissenseCount, firstImpactLikelihood, dndsLikelihood.missense());

        double substituteBoth = firstImpact == DriverImpact.MISSENSE || secondImpact == DriverImpact.MISSENSE
                ? multiHit
                : multiHit(nonBiallelicMissenseCount, nonBiallelicMissenseCount, dndsLikelihood.missense(), dndsLikelihood.missense());

        double combinedResult = Math.max(Math.max(substituteFirst, substituteSecond), substituteBoth);

        return builder.driverLikelihood(combinedResult).build();
    }

    private static double multiHit(int firstVariantTypeCount, int secondVariantTypeCount,
            final DndsDriverImpactLikelihood firstLikelihood, final DndsDriverImpactLikelihood secondLikelihood)
    {
        return probabilityDriverVariant(firstVariantTypeCount, secondVariantTypeCount, firstLikelihood, secondLikelihood);
    }

    private static double singleHit(int sampleCount, final DndsDriverImpactLikelihood likelihood)
    {
        return probabilityDriverVariant(sampleCount, likelihood);
    }

    private static int variantCount(
            boolean useBiallelic, final SomaticVariant variant,
            final Map<VariantType, Integer> standard, final Map<VariantType,Integer> biallelic)
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
