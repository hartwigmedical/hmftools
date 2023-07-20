package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.healthchecker.result.QCValue;

import org.jetbrains.annotations.NotNull;

public interface HealthChecker
{
    @NotNull
    List<QCValue> run() throws IOException;
}
