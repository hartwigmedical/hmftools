package com.hartwig.healthchecker.common.report;

import java.util.Optional;

import com.hartwig.healthchecker.common.exception.GenerateReportException;
import com.hartwig.healthchecker.common.io.dir.RunContext;
import com.hartwig.healthchecker.common.result.BaseResult;

import org.jetbrains.annotations.NotNull;

public interface Report {

    void addResult(@NotNull BaseResult result);

    @NotNull
    Optional<String> generateReport(@NotNull RunContext runContext, @NotNull String outputPath)
            throws GenerateReportException;
}
