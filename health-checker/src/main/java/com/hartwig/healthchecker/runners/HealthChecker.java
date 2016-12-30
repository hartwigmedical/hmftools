package com.hartwig.healthchecker.runners;

import com.hartwig.healthchecker.io.dir.RunContext;
import com.hartwig.healthchecker.result.BaseResult;

import org.jetbrains.annotations.NotNull;

public interface HealthChecker {

    @NotNull
    CheckType checkType();

    @NotNull
    BaseResult run(@NotNull RunContext runContext);
}
