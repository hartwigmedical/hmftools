package com.hartwig.hmftools.healthchecker.result;

import com.hartwig.hmftools.healthchecker.runners.CheckType;

import org.jetbrains.annotations.NotNull;

abstract class AbstractResult implements BaseResult {

    @NotNull
    private final CheckType checkType;

    AbstractResult(@NotNull final CheckType checkType) {
        this.checkType = checkType;
    }

    @Override
    @NotNull
    public CheckType getCheckType() {
        return checkType;
    }
}
