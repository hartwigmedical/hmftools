package com.hartwig.hmftools.healthchecker.result;

import com.hartwig.hmftools.healthchecker.runners.CheckType;

import org.jetbrains.annotations.NotNull;

public class NoResult extends AbstractResult {

    public NoResult(@NotNull final CheckType checkType) {
        super(checkType);
    }
}
