package com.hartwig.hmftools.protect.actionability;

import org.jetbrains.annotations.NotNull;

public interface Actionable {

    @NotNull
    String source();

    @NotNull
    String level();

    @NotNull
    String response();

    @NotNull
    String reference();

    @NotNull
    String drug();

    @NotNull
    String drugsType();

    @NotNull
    String cancerType();
}
