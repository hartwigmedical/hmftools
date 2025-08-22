package com.hartwig.hmftools.purple.germline;

import static java.util.stream.Collectors.toList;
import static java.util.stream.Collectors.toMap;
import static java.util.stream.Collectors.toSet;

import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.getWorstReportableCodingEffect;
import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.hasTranscriptCodingEffect;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.DriverImpact;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.CodingEffect;

import org.apache.logging.log4j.util.Strings;

public class GermlineDrivers
{
    private final Map<String, DriverCategory> mDriverCatalogMap;

    public GermlineDrivers(final List<DriverGene> driverGenes)
    {
        mDriverCatalogMap = driverGenes.stream()
                .filter(DriverGene::reportGermline).collect(toMap(DriverGene::gene, DriverGene::likelihoodType));
    }

    public List<DriverCatalog> findDrivers(final List<GermlineVariant> variants, final Map<String,List<GeneCopyNumber>> geneCopyNumberMap)
    {
        final Set<String> genes = variants.stream().map(GermlineVariant::gene).collect(toSet());
        final List<DriverCatalog> driverCatalog = Lists.newArrayList();

        for(String gene : genes)
        {
            DriverCategory category = mDriverCatalogMap.get(gene);
            List<GermlineVariant> geneVariants = variants.stream().filter(x -> x.gene().equals(gene)).collect(toList());

            if(category == null)
                continue;

            List<GeneCopyNumber> geneCopyNumbers = geneCopyNumberMap.get(gene);

            if(geneCopyNumbers == null)
                continue;

            for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
            {
                if(geneCopyNumbers.size() == 1)
                {
                    driverCatalog.add(germlineDriver(category, gene, geneVariants, geneCopyNumber));
                }
                else
                {
                    if(geneVariants.stream().anyMatch(x -> hasTranscriptCodingEffect(x.variantImpact(), x.type(), geneCopyNumber.TransName)))
                    {
                        driverCatalog.add(germlineDriver(category, gene, geneVariants, geneCopyNumber));
                    }
                }
            }
        }

        return driverCatalog;
    }

    private static DriverCatalog germlineDriver(
            final DriverCategory category, final String gene,
            final List<GermlineVariant> geneVariants, final GeneCopyNumber geneCopyNumber)
    {
        Map<DriverImpact,Integer> variantCounts = Maps.newHashMap();

        for(GermlineVariant variant : geneVariants)
        {
            CodingEffect codingEffect = getWorstReportableCodingEffect(variant.variantImpact());
            DriverImpact driverImpact = DriverImpact.select(variant.type(), codingEffect);
            Integer count = variantCounts.get(driverImpact);
            variantCounts.put(driverImpact, count != null ? count + 1 : 1);
        }

        int missenseVariants = variantCounts.getOrDefault(DriverImpact.MISSENSE, 0);
        int nonsenseVariants = variantCounts.getOrDefault(DriverImpact.NONSENSE, 0);
        int spliceVariants = variantCounts.getOrDefault(DriverImpact.SPLICE, 0);
        int inframeVariants = variantCounts.getOrDefault(DriverImpact.INFRAME, 0);
        int frameshiftVariants = variantCounts.getOrDefault(DriverImpact.FRAMESHIFT, 0);

        final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .chromosome(geneVariants.get(0).chromosome())
                .chromosomeBand(geneCopyNumber == null ? Strings.EMPTY : geneCopyNumber.ChromosomeBand)
                .gene(gene)
                .transcript(geneCopyNumber != null ? geneCopyNumber.TransName : "")
                .isCanonical(geneCopyNumber != null ? geneCopyNumber.IsCanonical : true)
                .driver(DriverType.GERMLINE_MUTATION)
                .category(category)
                .driverLikelihood(1)
                .missense(missenseVariants)
                .nonsense(nonsenseVariants)
                .splice(spliceVariants)
                .inframe(inframeVariants)
                .frameshift(frameshiftVariants)
                .biallelic(geneVariants.stream().anyMatch(GermlineVariant::biallelic))
                .minCopyNumber(geneCopyNumber == null ? 0 : geneCopyNumber.minCopyNumber())
                .maxCopyNumber(geneCopyNumber == null ? 0 : geneCopyNumber.maxCopyNumber())
                .likelihoodMethod(LikelihoodMethod.GERMLINE);

        return builder.build();
    }
}
