package com.hartwig.hmftools.purple.drivers;

import static java.util.stream.Collectors.counting;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toList;
import static java.util.stream.Collectors.toMap;
import static java.util.stream.Collectors.toSet;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverImpact;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
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

    public List<DriverCatalog> drivers(final List<VariantContext> variants, final Map<String,List<GeneCopyNumber>> geneCopyNumberMap)
    {
        final List<VariantContextDecorator> decoratedVariants = variants.stream().map(VariantContextDecorator::new).collect(toList());
        final Set<String> genes = decoratedVariants.stream().map(VariantContextDecorator::gene).collect(toSet());
        final List<DriverCatalog> driverCatalog = Lists.newArrayList();

        for(String gene : genes)
        {
            // TODO: decide how to handle multiple transcripts per gene
            GeneCopyNumber geneCopyNumber = geneCopyNumberMap.get(gene).get(0);

            DriverCategory category = mDriverCatalogMap.get(gene);
            List<VariantContextDecorator> geneVariants = decoratedVariants.stream().filter(x -> x.gene().equals(gene)).collect(toList());

            if(category != null)
                driverCatalog.add(germlineDriver(category, gene, geneVariants, geneCopyNumber));
        }

        return new ArrayList<>(driverCatalog);
    }

    static DriverCatalog germlineDriver(DriverCategory category, final String gene,
            final List<VariantContextDecorator> geneVariants, @Nullable GeneCopyNumber geneCopyNumber)
    {
        final Map<DriverImpact,Integer> variantCounts =
                geneVariants.stream().collect(groupingBy(VariantContextDecorator::impact, counting()))
                        .entrySet().stream().collect(Collectors.toMap(entry -> entry.getKey(), entry -> entry.getValue().intValue()));

        int missenseVariants = variantCounts.getOrDefault(DriverImpact.MISSENSE, 0);
        int nonsenseVariants = variantCounts.getOrDefault(DriverImpact.NONSENSE, 0);
        int spliceVariants = variantCounts.getOrDefault(DriverImpact.SPLICE, 0);
        int inframeVariants = variantCounts.getOrDefault(DriverImpact.INFRAME, 0);
        int frameshiftVariants = variantCounts.getOrDefault(DriverImpact.FRAMESHIFT, 0);

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
