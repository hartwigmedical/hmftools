package com.hartwig.healthchecker.result;

import com.hartwig.healthchecker.checkers.CheckType;

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
