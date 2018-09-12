package com.hartwig.hmftools.common.dnds;

import java.io.IOException;
import java.nio.charset.Charset;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;

public class DndsDriverLikelihoodSupplier {

    @NotNull
    @VisibleForTesting
    static Map<String, DndsDriverLikelihood> tsgLikelihood() throws IOException {
        return DndsDriverLikelihoodFile.fromLines(Resources.readLines(Resources.getResource("dnds/DndsDriverLikelihoodTsg.tsv"), Charset.defaultCharset()));
    }

    @NotNull
    @VisibleForTesting
    static Map<String, DndsDriverLikelihood> oncoLikelihood() throws IOException {
        return DndsDriverLikelihoodFile.fromLines(Resources.readLines(Resources.getResource("dnds/DndsDriverLikelihoodOnco.tsv"), Charset.defaultCharset()));
    }
}
