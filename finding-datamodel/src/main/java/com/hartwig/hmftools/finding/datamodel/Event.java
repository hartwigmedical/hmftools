package com.hartwig.hmftools.finding.datamodel;

import jakarta.validation.constraints.NotNull;

public interface Event extends Finding {

    @NotNull
    String event();
}
