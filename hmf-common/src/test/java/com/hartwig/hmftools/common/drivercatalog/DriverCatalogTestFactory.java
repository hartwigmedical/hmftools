package com.hartwig.hmftools.common.drivercatalog;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class DriverCatalogTestFactory {

    private DriverCatalogTestFactory() {
    }

    @NotNull
    public static DriverCatalog createNonCanonicalSomaticMutationEntryForGene(@NotNull String gene, double likelihood,
            @NotNull String transcript, @NotNull DriverCategory category) {
        return create(gene, likelihood, DriverType.MUTATION, transcript, false, category);
    }

    @NotNull
    public static DriverCatalog createCanonicalSomaticMutationEntryForGene(@NotNull String gene, double likelihood,
            @NotNull String transcript, @NotNull DriverCategory category) {
        return create(gene, likelihood, DriverType.MUTATION, transcript, true, category);
    }

    @NotNull
    public static DriverCatalog createCanonicalGermlineMutationEntryForGene(@NotNull String gene, double likelihood,
            @NotNull String transcript, @NotNull DriverCategory category) {
        return create(gene, likelihood, DriverType.GERMLINE_MUTATION, transcript, true, category);
    }

    @NotNull
    private static DriverCatalog create(@NotNull String gene, double likelihood, @NotNull DriverType type, @NotNull String transcript,
            boolean isCanonical, @NotNull DriverCategory category) {
        return ImmutableDriverCatalog.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(gene)
                .transcript(transcript)
                .isCanonical(isCanonical)
                .driver(type)
                .category(category)
                .likelihoodMethod(LikelihoodMethod.DNDS)
                .driverLikelihood(likelihood)
                .missense(0)
                .nonsense(0)
                .splice(0)
                .inframe(0)
                .frameshift(0)
                .biallelic(false)
                .minCopyNumber(0)
                .maxCopyNumber(0)
                .build();
    }
}