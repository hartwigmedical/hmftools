package com.hartwig.hmftools.common.drivercatalog;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.dnds.DndsDriverLikelihood;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.jetbrains.annotations.NotNull;

public class OncoDrivers {

    private final Map<String, DndsDriverLikelihood> refCdsCv;

    public OncoDrivers(@NotNull final Map<String, DndsDriverLikelihood> refCdsCv) {
        this.refCdsCv = refCdsCv;
    }

    public Set<DriverCatalog> doSTuff(@NotNull final List<DndsDriverLikelihood> likelihoods,
            @NotNull final List<EnrichedSomaticVariant> variants) {

        final Set<DriverCatalog> driverCatalog = Sets.newHashSet();
        final Map<String, DndsDriverLikelihood> likelihoodsByGene =
                likelihoods.stream().collect(Collectors.toMap(DndsDriverLikelihood::gene, x -> x));

        final Map<String, List<EnrichedSomaticVariant>> missenseVariantsByGene = variants.stream()
                .filter(x -> likelihoodsByGene.keySet().contains(x.gene()))
                .filter(x -> x.canonicalCodingEffect().equals(CodingEffect.MISSENSE))
                .collect(Collectors.groupingBy(SomaticVariant::gene));

        final Map<VariantType, Long> variantTypeCounts =
                variants.stream().collect(Collectors.groupingBy(SomaticVariant::type, Collectors.counting()));
        long sampleSNVCount = variantTypeCounts.getOrDefault(VariantType.SNP, 0L);

        for (String gene : missenseVariantsByGene.keySet()) {
            final ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder().gene(gene).likelihood(1);

            final DndsDriverLikelihood likelihood = likelihoodsByGene.get(gene);
            final double missenseDriverLikelihood = missenseProbabilityDriverVariant(sampleSNVCount, likelihood);

            final List<EnrichedSomaticVariant> geneVariants = missenseVariantsByGene.get(gene);
            if (geneVariants.stream().anyMatch(SomaticVariant::hotspot)) {
                driverCatalog.add(builder.driver(DriverType.HOTSPOT).build());
            }

            final List<EnrichedSomaticVariant> inframe =
                    geneVariants.stream().filter(x -> x.type() == VariantType.INDEL && x.repeatCount() < 8).collect(Collectors.toList());
            if (!inframe.isEmpty()) {
                driverCatalog.add(builder.driver(DriverType.INFRAME).build());
            }

            final List<EnrichedSomaticVariant> pointMutations =
                    geneVariants.stream().filter(x -> x.type() != VariantType.INDEL).collect(Collectors.toList());
            if (Doubles.greaterThan(missenseDriverLikelihood, 0) && pointMutations.stream().anyMatch(x -> !x.hotspot())) {
                driverCatalog.add(builder.driver(geneVariants.size() > 1 ? DriverType.MULTI_HIT : DriverType.SINGLE_HIT)
                        .likelihood(missenseDriverLikelihood)
                        .build());
            }
        }

        return driverCatalog;
    }

    public static double missenseProbabilityDriverVariant(long sampleSNVCount, @NotNull final DndsDriverLikelihood likelihood) {

        double lambda = 1.0 * sampleSNVCount / likelihood.missenseProbabilityVariantNonDriverFactor();
        PoissonDistribution poissonDistribution = new PoissonDistribution(lambda);

        double pDriver = likelihood.missenseProbabilityDriver();
        double pVariantNondriver = 1 - poissonDistribution.cumulativeProbability(0);

        return pDriver / (pDriver + pVariantNondriver * (1 - pDriver));
    }

}
