package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.drivercatalog.DriverCatalogFactory.probabilityDriverVariant;

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
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class OncoDrivers
{
    static final int MAX_REPEAT_COUNT = 7;

    private final ReportablePredicate mReportablePredicate;
    private final Map<String, DndsDriverGeneLikelihood> mLikelihoodsByGene;

    public OncoDrivers(@NotNull DriverGenePanel genePanel)
    {
        mLikelihoodsByGene = genePanel.oncoLikelihood();
        mReportablePredicate = new ReportablePredicate(DriverCategory.ONCO, genePanel);
    }

    @NotNull
    List<DriverCatalog> drivers(@NotNull final List<SomaticVariant> variants, @NotNull final List<GeneCopyNumber> geneCopyNumberList,
            @NotNull final Map<VariantType, Long> variantTypeCounts)
    {
        final List<DriverCatalog> driverCatalog = Lists.newArrayList();

        final Map<String, GeneCopyNumber> geneCopyNumbers =
                geneCopyNumberList.stream().collect(Collectors.toMap(TranscriptRegion::geneName, x -> x));
        long sampleSNVCount = variantTypeCounts.getOrDefault(VariantType.SNP, 0L);
        long sampleINDELCount = variantTypeCounts.getOrDefault(VariantType.INDEL, 0L);

        final Map<String, List<SomaticVariant>> codingVariants = oncogenicVariantsByGene(variants);

        for(String gene : codingVariants.keySet())
        {
            final DndsDriverGeneLikelihood geneMissenseLikelihood = mLikelihoodsByGene.get(gene);

            final List<SomaticVariant> geneVariants = codingVariants.get(gene);

            driverCatalog.add(geneDriver(
                    sampleSNVCount, sampleINDELCount, gene, geneMissenseLikelihood, geneVariants, geneCopyNumbers.get(gene)));
        }

        return driverCatalog;
    }

    @NotNull
    private Map<String, List<SomaticVariant>> oncogenicVariantsByGene(@NotNull final List<SomaticVariant> variants)
    {
        return variants.stream().filter(mReportablePredicate).collect(Collectors.groupingBy(SomaticVariant::gene));
    }

    @NotNull
    static DriverCatalog geneDriver(long sampleSNVCount, long sampleIndelCount, @NotNull final String gene,
            @Nullable final DndsDriverGeneLikelihood geneLikelihood, @NotNull final List<SomaticVariant> geneVariants,
            @Nullable GeneCopyNumber geneCopyNumber)
    {
        final Map<DriverImpact, Long> variantCounts = DriverCatalogFactory.driverImpactCount(geneVariants);
        long missenseVariants = variantCounts.getOrDefault(DriverImpact.MISSENSE, 0L);
        long nonsenseVariants = variantCounts.getOrDefault(DriverImpact.NONSENSE, 0L);
        long spliceVariants = variantCounts.getOrDefault(DriverImpact.SPLICE, 0L);
        long inframeVariants = variantCounts.getOrDefault(DriverImpact.INFRAME, 0L);
        long frameshiftVariants = variantCounts.getOrDefault(DriverImpact.FRAMESHIFT, 0L);

        final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .chromosome(geneVariants.get(0).chromosome())
                .chromosomeBand(geneCopyNumber == null ? Strings.EMPTY : geneCopyNumber.chromosomeBand())
                .gene(gene)
                .driver(DriverType.MUTATION)
                .category(DriverCategory.ONCO)
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

        if(geneVariants.stream().anyMatch(OncoDrivers::isKnownInframeIndel))
        {
            return builder.likelihoodMethod(LikelihoodMethod.INFRAME).build();
        }

        double driverLikelihood = 0;

        if(geneLikelihood != null)
        {
            for(SomaticVariant variant : geneVariants)
            {
                final DriverImpact impact = DriverImpact.select(variant);

                final DndsDriverImpactLikelihood likelihood = geneLikelihood.select(impact);

                final long sampleVariantCount =
                        impact == DriverImpact.FRAMESHIFT || impact == DriverImpact.INFRAME ? sampleIndelCount : sampleSNVCount;

                driverLikelihood = Math.max(driverLikelihood, DriverCatalogFactory.probabilityDriverVariant(sampleVariantCount, likelihood));
            }
        }

        return builder.driverLikelihood(driverLikelihood).build();
    }

    private static boolean isKnownInframeIndel(@NotNull final SomaticVariant variant)
    {
        return variant.type() == VariantType.INDEL && variant.canonicalCodingEffect() == CodingEffect.MISSENSE
                && variant.repeatCount() <= MAX_REPEAT_COUNT;
    }
}
