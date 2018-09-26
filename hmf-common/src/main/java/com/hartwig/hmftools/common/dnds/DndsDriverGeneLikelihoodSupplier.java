package com.hartwig.hmftools.common.dnds;

import java.io.InputStream;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

public final class DndsDriverGeneLikelihoodSupplier {

    private DndsDriverGeneLikelihoodSupplier() {
    }

    @NotNull
    @VisibleForTesting
    public static Map<String, DndsDriverGeneLikelihood> tsgLikelihood() {
        final InputStream inputStream = DndsDriverGeneLikelihoodSupplier.class.getResourceAsStream("/dnds/DndsDriverLikelihoodTsg.tsv");
        return DndsDriverGeneLikelihoodFile.fromInputStream(inputStream);
    }

    @NotNull
    @VisibleForTesting
    public static Map<String, DndsDriverGeneLikelihood> oncoLikelihood() {
        final InputStream inputStream = DndsDriverGeneLikelihoodSupplier.class.getResourceAsStream("/dnds/DndsDriverLikelihoodOnco.tsv");
        return DndsDriverGeneLikelihoodFile.fromInputStream(inputStream);
    }
}
