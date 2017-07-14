package com.hartwig.hmftools.common.ecrf.reader;

import org.jetbrains.annotations.NotNull;

public interface OIDObject {

    @NotNull
    String OID();

    @NotNull
    String name();
}
