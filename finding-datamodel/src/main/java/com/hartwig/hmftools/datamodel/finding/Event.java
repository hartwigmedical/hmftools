package com.hartwig.hmftools.datamodel.finding;

import jakarta.validation.constraints.NotNull;

public interface Event extends Finding {

    @NotNull
    String event();
}
