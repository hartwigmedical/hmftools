package com.hartwig.hmftools.common.drivercatalog;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.region.TranscriptRegion;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class GermlineDrivers {

    private final Map<String, DriverCategory> driverCatalogMap;

    public GermlineDrivers(final List<DriverGene> driverGenes) {
        driverCatalogMap = driverGenes.stream()
                .filter(DriverGene::reportGermline)
                .collect(Collectors.toMap(DriverGene::gene, DriverGene::likelihoodType));
    }

    @NotNull
    public List<DriverCatalog> drivers(@NotNull final List<VariantContext> variants,
            @NotNull final List<GeneCopyNumber> geneCopyNumberList) {
        final Map<String, GeneCopyNumber> geneCopyNumbers =
                geneCopyNumberList.stream().collect(Collectors.toMap(TranscriptRegion::gene, x -> x));
        final Set<DriverCatalog> driverCatalog = Sets.newHashSet();

        for (VariantContext variant : variants) {
            final VariantContextDecorator decorator = new VariantContextDecorator(variant);
            final String gene = decorator.gene();
            if (decorator.reported()) {
                DriverCategory category = driverCatalogMap.get(gene);
                GeneCopyNumber geneCopyNumber = geneCopyNumbers.get(gene);
                if (category != null && geneCopyNumber != null) {
                    driverCatalog.add(cnaDriver(category, geneCopyNumber, decorator.biallelic()));
                }
            }
        }
        return new ArrayList<>(driverCatalog);
    }

    @NotNull
    private DriverCatalog cnaDriver(DriverCategory category, GeneCopyNumber x, boolean isBiallelic) {
        return ImmutableDriverCatalog.builder()
                .chromosome(x.chromosome())
                .chromosomeBand(x.chromosomeBand())
                .gene(x.gene())
                .missense(0)
                .nonsense(0)
                .inframe(0)
                .frameshift(0)
                .splice(0)
                .dndsLikelihood(0)
                .driverLikelihood(1)
                .driver(DriverType.MUTATION)
                .likelihoodMethod(LikelihoodMethod.GERMLINE)
                .category(category)
                .biallelic(isBiallelic)
                .minCopyNumber(x.minCopyNumber())
                .maxCopyNumber(x.maxCopyNumber())
                .build();
    }

}
