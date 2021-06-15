package com.hartwig.hmftools.common.drivercatalog;

import static java.util.stream.Collectors.counting;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toList;
import static java.util.stream.Collectors.toMap;
import static java.util.stream.Collectors.toSet;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.region.TranscriptRegion;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

public class GermlineDrivers
{
    @NotNull
    private final Map<String, DriverCategory> mDriverCatalogMap;

    public GermlineDrivers(@NotNull final List<DriverGene> driverGenes)
    {
        mDriverCatalogMap = driverGenes.stream()
                .filter(DriverGene::reportGermline).collect(toMap(DriverGene::gene, DriverGene::likelihoodType));
    }

    @NotNull
    public List<DriverCatalog> drivers(@NotNull final List<VariantContext> variants,
            @NotNull final List<GeneCopyNumber> geneCopyNumberList)
    {
        final List<VariantContextDecorator> decoratedVariants = variants.stream().map(VariantContextDecorator::new).collect(toList());
        final Set<String> genes = decoratedVariants.stream().map(VariantContextDecorator::gene).collect(toSet());
        final Map<String, GeneCopyNumber> geneCopyNumbers = geneCopyNumberList.stream().collect(toMap(TranscriptRegion::gene, x -> x));
        final List<DriverCatalog> driverCatalog = Lists.newArrayList();
        for(String gene : genes)
        {
            GeneCopyNumber geneCopyNumber = geneCopyNumbers.get(gene);
            DriverCategory category = mDriverCatalogMap.get(gene);
            List<VariantContextDecorator> geneVariants = decoratedVariants.stream().filter(x -> x.gene().equals(gene)).collect(toList());
            if(category != null)
            {
                driverCatalog.add(germlineDriver(category, gene, geneVariants, geneCopyNumber));
            }
        }

        return new ArrayList<>(driverCatalog);
    }

    @NotNull
    static DriverCatalog germlineDriver(DriverCategory category, @NotNull final String gene,
            @NotNull final List<VariantContextDecorator> geneVariants, @Nullable GeneCopyNumber geneCopyNumber)
    {
        final Map<DriverImpact, Long> variantCounts =
                geneVariants.stream().collect(groupingBy(VariantContextDecorator::impact, counting()));
        long missenseVariants = variantCounts.getOrDefault(DriverImpact.MISSENSE, 0L);
        long nonsenseVariants = variantCounts.getOrDefault(DriverImpact.NONSENSE, 0L);
        long spliceVariants = variantCounts.getOrDefault(DriverImpact.SPLICE, 0L);
        long inframeVariants = variantCounts.getOrDefault(DriverImpact.INFRAME, 0L);
        long frameshiftVariants = variantCounts.getOrDefault(DriverImpact.FRAMESHIFT, 0L);

        final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .chromosome(geneVariants.get(0).chromosome())
                .chromosomeBand(geneCopyNumber == null ? Strings.EMPTY : geneCopyNumber.chromosomeBand())
                .gene(gene)
                .driver(DriverType.GERMLINE)
                .category(category)
                .driverLikelihood(1)
                .missense(missenseVariants)
                .nonsense(nonsenseVariants)
                .splice(spliceVariants)
                .inframe(inframeVariants)
                .frameshift(frameshiftVariants)
                .biallelic(geneVariants.stream().anyMatch(VariantContextDecorator::biallelic))
                .minCopyNumber(geneCopyNumber == null ? 0 : geneCopyNumber.minCopyNumber())
                .maxCopyNumber(geneCopyNumber == null ? 0 : geneCopyNumber.maxCopyNumber())
                .likelihoodMethod(LikelihoodMethod.GERMLINE);

        return builder.build();
    }
}
