package com.hartwig.hmftools.common.dnds;

import java.io.IOException;
import java.nio.charset.Charset;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;

public class DndsDriverGeneLikelihoodSupplier {

    @NotNull
    @VisibleForTesting
    public static Map<String, DndsDriverGeneLikelihood> tsgLikelihood() throws IOException {
        return DndsDriverGeneLikelihoodFile.fromLines(Resources.readLines(Resources.getResource("dnds/DndsDriverLikelihoodTsg.tsv"),
                Charset.defaultCharset()));
    }

    @NotNull
    @VisibleForTesting
    public static Map<String, DndsDriverGeneLikelihood> oncoLikelihood() throws IOException {
        return DndsDriverGeneLikelihoodFile.fromLines(Resources.readLines(Resources.getResource("dnds/DndsDriverLikelihoodOnco.tsv"),
                Charset.defaultCharset()));
    }
}
