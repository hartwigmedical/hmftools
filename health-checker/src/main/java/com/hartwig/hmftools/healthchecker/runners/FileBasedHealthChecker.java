package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;
import java.util.function.Function;

import com.hartwig.hmftools.healthchecker.result.QCValue;

import org.jetbrains.annotations.NotNull;

public class FileBasedHealthChecker<T> implements HealthChecker {

    @NotNull
    private final T fileToQC;
    @NotNull
    private final Function<T, List<QCValue>> qcValueFunction;

    public FileBasedHealthChecker(@NotNull T fileToQC, @NotNull Function<T, List<QCValue>> qcValueFunction) {
        this.fileToQC = fileToQC;
        this.qcValueFunction = qcValueFunction;
    }

    @NotNull
    @Override
    public List<QCValue> run() throws IOException {
        return qcValueFunction.apply(fileToQC);
    }
}
