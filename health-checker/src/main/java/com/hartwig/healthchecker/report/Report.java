package com.hartwig.healthchecker.report;

import java.util.Optional;

import com.hartwig.healthchecker.exception.GenerateReportException;
import com.hartwig.healthchecker.io.dir.RunContext;
import com.hartwig.healthchecker.result.BaseResult;

import org.jetbrains.annotations.NotNull;

public interface Report {

    void addResult(@NotNull BaseResult result);

    @NotNull
    Optional<String> generateReport(@NotNull RunContext runContext, @NotNull String outputPath)
            throws GenerateReportException;
}
