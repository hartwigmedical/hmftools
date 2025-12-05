package com.hartwig.hmftools.common.driver;

import com.hartwig.hmftools.common.purple.ReportedStatus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class DriverCatalogTestFactory {

    private DriverCatalogTestFactory() {
    }

    public static DriverCatalog createCanonicalSomaticMutationEntryForGene(
            final String gene, double likelihood,
            final String transcript, final DriverCategory category)
    {
        return builder().gene(gene)
                .transcript(transcript)
                .isCanonical(true)
                .driverLikelihood(likelihood)
                .driver(DriverType.MUTATION)
                .category(category)
                .build();
    }

    @NotNull
    public static ImmutableDriverCatalog.Builder builder()
    {
        return ImmutableDriverCatalog.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .isCanonical(false)
                .driver(DriverType.MUTATION)
                .category(DriverCategory.TSG)
                .reportedStatus(ReportedStatus.REPORTED)
                .likelihoodMethod(LikelihoodMethod.DNDS)
                .driverLikelihood(0D)
                .missense(0)
                .nonsense(0)
                .splice(0)
                .inframe(0)
                .frameshift(0)
                .biallelic(false)
                .minCopyNumber(0)
                .maxCopyNumber(0);
    }
}