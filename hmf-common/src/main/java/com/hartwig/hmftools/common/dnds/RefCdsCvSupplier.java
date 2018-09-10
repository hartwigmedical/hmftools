package com.hartwig.hmftools.common.dnds;

import java.io.IOException;
import java.nio.charset.Charset;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;

public class RefCdsCvSupplier {

    @NotNull
    @VisibleForTesting
    static Map<String, RefCdsCv> refCdsCv() throws IOException {
        return RefCdsCvFile.fromLines(Resources.readLines(Resources.getResource("dnds/HmfRefCdsCv.csv"), Charset.defaultCharset()));
    }
}
