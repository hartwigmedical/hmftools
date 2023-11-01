package com.hartwig.hmftools.healthchecker.runners;

import java.util.List;
import java.util.function.Function;

import com.hartwig.hmftools.healthchecker.result.QCValue;

public class FileBasedHealthChecker<T> implements HealthChecker
{
    private final T fileToQC;
    private final Function<T, List<QCValue>> qcValueFunction;

    public FileBasedHealthChecker(final T fileToQC, final Function<T, List<QCValue>> qcValueFunction)
    {
        this.fileToQC = fileToQC;
        this.qcValueFunction = qcValueFunction;
    }

    @Override
    public List<QCValue> run()
    {
        return qcValueFunction.apply(fileToQC);
    }
}
