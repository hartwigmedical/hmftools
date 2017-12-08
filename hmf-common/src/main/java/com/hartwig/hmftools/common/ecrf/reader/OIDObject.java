package com.hartwig.hmftools.common.ecrf.reader;

import org.jetbrains.annotations.NotNull;

interface OIDObject {

    @NotNull
    String oid();

    @NotNull
    String name();
}
