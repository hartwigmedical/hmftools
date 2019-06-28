package com.hartwig.hmftools.common.dnds;

import static com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodFile.fromMultiImpactInputStream;
import static com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodFile.fromSingleImpactInputStream;

import java.io.InputStream;
import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public final class DndsDriverGeneLikelihoodSupplier {

    private DndsDriverGeneLikelihoodSupplier() {
    }

    @NotNull
    public static Map<String, DndsDriverGeneLikelihood> tsgLikelihood() {
        final Map<String, DndsDriverGeneLikelihood> result = Maps.newHashMap();

        final Map<String, ModifiableDndsDriverGeneLikelihood> standard =
                fromMultiImpactInputStream(resource("/dnds/DndsDriverLikelihoodTsg.tsv"));
        final Map<String, DndsDriverImpactLikelihood> biallelic =
                fromSingleImpactInputStream(resource("/dnds/DndsDriverLikelihoodTsgBiallelic.tsv"));
        final Map<String, DndsDriverImpactLikelihood> nonBiallelic =
                fromSingleImpactInputStream(resource("/dnds/DndsDriverLikelihoodTsgNonBiallelic.tsv"));

        for (ModifiableDndsDriverGeneLikelihood likelihood : standard.values()) {
            final String gene = likelihood.gene();
            if (biallelic.keySet().contains(gene)) {
                likelihood.setMissenseBiallelic(biallelic.get(gene));
            }

            if (nonBiallelic.keySet().contains(gene)) {
                likelihood.setMissenseNonBiallelic(nonBiallelic.get(gene));
            }

            result.put(gene, likelihood);
        }

        return result;
    }

    @NotNull
    public static Map<String, DndsDriverImpactLikelihood> oncoLikelihood() {
        return fromSingleImpactInputStream(resource("/dnds/DndsDriverLikelihoodOnco.tsv"));
    }

    @NotNull
    private static InputStream resource(@NotNull final String resource) {
        return DndsDriverGeneLikelihoodSupplier.class.getResourceAsStream(resource);
    }

}
