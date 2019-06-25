package com.hartwig.hmftools.common.dnds;

import java.io.InputStream;
import java.util.Map;

import org.jetbrains.annotations.NotNull;

public final class DndsDriverGeneLikelihoodSupplier {

    private DndsDriverGeneLikelihoodSupplier() {
    }

    @NotNull
    public static Map<String, DndsDriverGeneLikelihood> tsgLikelihood() {
        final InputStream inputStream = DndsDriverGeneLikelihoodSupplier.class.getResourceAsStream("/dnds/DndsDriverLikelihoodTsg.tsv");
        return DndsDriverGeneLikelihoodFile.fromTsgInputStream(inputStream);
    }

    @NotNull
    public static Map<String, DndsDriverImpactLikelihood> oncoLikelihood() {
        final InputStream inputStream = DndsDriverGeneLikelihoodSupplier.class.getResourceAsStream("/dnds/DndsDriverLikelihoodOnco.tsv");
        return DndsDriverGeneLikelihoodFile.fromOncoInputStream(inputStream);
    }
}
