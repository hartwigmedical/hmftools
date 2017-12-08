package com.hartwig.hmftools.healthchecker.result;

import com.hartwig.hmftools.healthchecker.runners.CheckType;

import org.jetbrains.annotations.NotNull;

public interface BaseResult {

    @NotNull
    CheckType getCheckType();
}
