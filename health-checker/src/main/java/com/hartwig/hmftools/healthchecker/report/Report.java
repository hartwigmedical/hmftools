package com.hartwig.hmftools.healthchecker.report;

import java.io.IOException;
import java.util.Optional;

import com.hartwig.hmftools.healthchecker.result.BaseResult;

import org.jetbrains.annotations.NotNull;

public interface Report {

    void addResult(@NotNull BaseResult result);

    @NotNull
    Optional<String> generateReport(@NotNull String fileName) throws IOException;
}
