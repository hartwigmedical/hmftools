package com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel;

import org.jetbrains.annotations.NotNull;

public class EcrfResolveException extends Exception {

    public EcrfResolveException(@NotNull String message) {
        super(message);
    }
}
