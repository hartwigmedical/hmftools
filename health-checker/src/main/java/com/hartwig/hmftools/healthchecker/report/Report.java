package com.hartwig.hmftools.healthchecker.report;

import java.util.Optional;

import com.hartwig.hmftools.common.exception.GenerateReportException;
import com.hartwig.hmftools.healthchecker.context.RunContext;
import com.hartwig.hmftools.healthchecker.result.BaseResult;

import org.jetbrains.annotations.NotNull;

public interface Report {

    void addResult(@NotNull BaseResult result);

    @NotNull
    Optional<String> generateReport(@NotNull RunContext runContext, @NotNull String outputPath)
            throws GenerateReportException;
}
