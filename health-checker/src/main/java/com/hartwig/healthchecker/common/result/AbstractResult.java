package com.hartwig.healthchecker.common.result;

import com.hartwig.healthchecker.common.checks.CheckType;

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
