package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.purple.PurpleConstants.MAX_INDEL_DRIVER_REPEAT_COUNT;
import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.getWorstReportableCodingEffect;
import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.groupByImpact;
import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.hasTranscriptCodingEffect;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogFactory;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.DriverImpact;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.driver.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.driver.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanel;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.purple.DriverSourceData;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

public class OncoDrivers extends SomaticVariantDriverFinder
{
    public OncoDrivers(final DriverGenePanel genePanel)
    {
        super(genePanel, DriverCategory.ONCO);
    }
    
    public List<DriverCatalog> findDrivers(
            final Map<String,List<GeneCopyNumber>> geneCopyNumberMap, final Map<VariantType,Integer> variantTypeCounts,
            final Map<VariantType,Integer> variantTypeCountsBiallelic, final List<DriverSourceData> driverSourceData)
    {
        List<DriverCatalog> driverCatalog = Lists.newArrayList();

        int sampleSNVCount = variantTypeCounts.getOrDefault(VariantType.SNP, 0);
        int sampleINDELCount = variantTypeCounts.getOrDefault(VariantType.INDEL, 0);

        final Map<String,List<SomaticVariant>> codingVariants = mReportableVariants.stream()
                .collect(Collectors.groupingBy(SomaticVariant::gene));

        for(String gene : codingVariants.keySet())
        {
            DndsDriverGeneLikelihood dndsLikelihood = mLikelihoodsByGene.containsKey(gene) ?
                    mLikelihoodsByGene.get(gene) : NO_GENE_DNDS_LIKELIHOOD;

            List<SomaticVariant> geneVariants = codingVariants.get(gene);

            List<GeneCopyNumber> geneCopyNumbers = geneCopyNumberMap.get(gene);

            if(geneCopyNumbers == null)
                continue;

            for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
            {
                if(geneCopyNumbers.size() == 1
                || geneVariants.stream().anyMatch(x -> hasTranscriptCodingEffect(x.variantImpact(), x.type(), geneCopyNumber.TransName)))
                {
                    // confirm this variant has a reportable effect against the specific transcript or the gene itself
                    DriverCatalog driverRecord = createOncoDriver(
                            sampleSNVCount, sampleINDELCount, dndsLikelihood, geneVariants, geneCopyNumber);
                    driverCatalog.add(driverRecord);

                    driverSourceData.add(new DriverSourceData(driverRecord, geneVariants.get(0)));
                }
            }
        }

        return driverCatalog;
    }

    public static DriverCatalog createOncoDriver(
            int sampleSNVCount, int sampleIndelCount, final DndsDriverGeneLikelihood geneLikelihood,
            final List<SomaticVariant> geneVariants, final GeneCopyNumber geneCopyNumber)
    {
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
                CodingEffect codingEffect = getWorstReportableCodingEffect(variant.variantImpact());
                final DriverImpact impact = DriverImpact.select(variant.type(), codingEffect);

                final DndsDriverImpactLikelihood likelihood = geneLikelihood.select(impact);

                final int sampleVariantCount =
                        impact == DriverImpact.FRAMESHIFT || impact == DriverImpact.INFRAME ? sampleIndelCount : sampleSNVCount;

                driverLikelihood = Math.max(driverLikelihood, DriverCatalogFactory.probabilityDriverVariant(sampleVariantCount, likelihood));
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
