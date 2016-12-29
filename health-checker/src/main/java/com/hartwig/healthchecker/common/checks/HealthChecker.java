package com.hartwig.healthchecker.common.checks;

import com.hartwig.healthchecker.common.io.dir.RunContext;
import com.hartwig.healthchecker.common.result.BaseResult;

import org.jetbrains.annotations.NotNull;

public interface HealthChecker {

    @NotNull
    CheckType checkType();

    @NotNull
    BaseResult run(@NotNull RunContext runContext);
}
