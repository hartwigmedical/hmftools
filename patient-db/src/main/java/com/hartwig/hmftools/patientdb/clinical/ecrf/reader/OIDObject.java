package com.hartwig.hmftools.patientdb.clinical.ecrf.reader;

import org.jetbrains.annotations.NotNull;

interface OIDObject {

    @NotNull
    String oid();

    @NotNull
    String name();
}
