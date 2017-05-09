package com.hartwig.hmftools.healthchecker.runners;

import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.healthchecker.result.BaseResult;

import org.jetbrains.annotations.NotNull;

public interface HealthChecker {

    @NotNull
    CheckType checkType();

    @NotNull
    BaseResult run(@NotNull RunContext runContext);
}
