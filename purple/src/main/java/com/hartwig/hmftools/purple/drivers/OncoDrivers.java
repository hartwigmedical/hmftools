package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.purple.PurpleConstants.MAX_INDEL_DRIVER_REPEAT_COUNT;
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

public class OncoDrivers extends SomaticVariantDriverFinder
{
    public OncoDrivers(final DriverGeneResource driverGeneResource)
    {
        super(driverGeneResource, DriverCategory.ONCO);
    }

    public DriverCatalog createDriverCatalog(
            final List<SomaticVariant> geneVariants, final Map<VariantType,Integer> variantTypeCounts,
            final Map<VariantType,Integer> biallelicCounts, final GeneCopyNumber geneCopyNumber,
            final DndsDriverGeneLikelihood dndsLikelihood)
    {
        int sampleSnvCount = variantTypeCounts.getOrDefault(VariantType.SNP, 0);
        int sampleIndelCount = variantTypeCounts.getOrDefault(VariantType.INDEL, 0);

        Map<DriverImpact,Integer> variantCounts = groupByImpact(geneVariants);
        int missenseVariants = variantCounts.getOrDefault(DriverImpact.MISSENSE, 0);
        int nonsenseVariants = variantCounts.getOrDefault(DriverImpact.NONSENSE, 0);
        int spliceVariants = variantCounts.getOrDefault(DriverImpact.SPLICE, 0);
        int inframeVariants = variantCounts.getOrDefault(DriverImpact.INFRAME, 0);
        int frameshiftVariants = variantCounts.getOrDefault(DriverImpact.FRAMESHIFT, 0);

        final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .chromosome(geneVariants.get(0).chromosome())
                .chromosomeBand(geneCopyNumber.ChromosomeBand)
                .gene(geneCopyNumber.geneName())
                .transcript(geneCopyNumber.TransName)
                .isCanonical(geneCopyNumber.IsCanonical)
                .driver(DriverType.MUTATION)
                .category(DriverCategory.ONCO)
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

        if(geneVariants.stream().anyMatch(OncoDrivers::isKnownInframeIndel))
        {
            return builder.likelihoodMethod(LikelihoodMethod.INFRAME).build();
        }

        double driverLikelihood = 0;

        if(dndsLikelihood != null)
        {
            for(SomaticVariant variant : geneVariants)
            {
                CodingEffect codingEffect = getWorstReportableCodingEffect(variant.variantImpact());
                DriverImpact impact = DriverImpact.select(variant.type(), codingEffect);

                DndsDriverImpactLikelihood likelihood = dndsLikelihood.select(impact);

                int sampleVariantCount = (impact == DriverImpact.FRAMESHIFT || impact == DriverImpact.INFRAME)
                        ? sampleIndelCount : sampleSnvCount;

                driverLikelihood = Math.max(driverLikelihood, DndsCalculator.probabilityDriverVariant(sampleVariantCount, likelihood));
            }
        }

        return builder.driverLikelihood(driverLikelihood).build();
    }

    private static boolean isKnownInframeIndel(final SomaticVariant variant)
    {
        return variant.type() == VariantType.INDEL && variant.variantImpact().CanonicalCodingEffect == CodingEffect.MISSENSE
                && variant.decorator().repeatCount() <= MAX_INDEL_DRIVER_REPEAT_COUNT;
    }
}
