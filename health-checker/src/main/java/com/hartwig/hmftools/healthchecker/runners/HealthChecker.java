package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;

import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.healthchecker.result.BaseResult;

import org.jetbrains.annotations.NotNull;

public interface HealthChecker {

    @NotNull
    BaseResult run(@NotNull RunContext runContext) throws IOException;
}
