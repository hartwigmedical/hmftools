package com.hartwig.hmftools.common.drivercatalog.dnds;

import static com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihoodFile.fromMultiImpactInputStream;
import static com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihoodFile.fromSingleImpactInputStream;

import java.io.InputStream;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public final class DndsDriverGeneLikelihoodSupplier {

    private DndsDriverGeneLikelihoodSupplier() {
    }

    @NotNull
    public static Map<String, DndsDriverGeneLikelihood> tsgLikelihood(@NotNull final Set<String> tsGenes) {
        Map<String, DndsDriverGeneLikelihood> filtered = Maps.newHashMap();
        Map<String, DndsDriverGeneLikelihood> all = tsgLikelihood();

        for (String target : tsGenes) {
            if (!all.containsKey(target)) {
                throw new IllegalArgumentException(target + " is not a valid TSG driver gene.");
            }
            filtered.put(target, all.get(target));
        }

        return filtered;
    }

    @NotNull
    public static Map<String, DndsDriverGeneLikelihood> oncoLikelihood(@NotNull final Set<String> tsGenes) {
        Map<String, DndsDriverGeneLikelihood> filtered = Maps.newHashMap();
        Map<String, ModifiableDndsDriverGeneLikelihood> all = oncoLikelihood();

        for (String target : tsGenes) {
            if (!all.containsKey(target)) {
                throw new IllegalArgumentException(target + " is not a valid ONCO driver gene.");
            }
            filtered.put(target, all.get(target));
        }

        return filtered;
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
            if (biallelic.containsKey(gene)) {
                likelihood.setMissenseBiallelic(biallelic.get(gene));
            }

            if (nonBiallelic.containsKey(gene)) {
                likelihood.setMissenseNonBiallelic(nonBiallelic.get(gene));
            }

            result.put(gene, likelihood);
        }

        return result;
    }

    @NotNull
    public static Map<String, ModifiableDndsDriverGeneLikelihood> oncoLikelihood() {
        return fromMultiImpactInputStream(resource("/dnds/DndsDriverLikelihoodOnco.tsv"));
    }

    @NotNull
    private static InputStream resource(@NotNull final String resource) {
        return DndsDriverGeneLikelihoodSupplier.class.getResourceAsStream(resource);
    }
}
