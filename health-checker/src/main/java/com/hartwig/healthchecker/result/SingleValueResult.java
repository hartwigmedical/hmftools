package com.hartwig.healthchecker.result;

import com.hartwig.healthchecker.runners.CheckType;
import com.hartwig.healthchecker.runners.checks.HealthCheck;

import org.jetbrains.annotations.NotNull;

public class SingleValueResult extends AbstractResult {

    private static final long serialVersionUID = -5744830259786248569L;

    @NotNull
    private final HealthCheck check;

    public SingleValueResult(@NotNull final CheckType checkType, @NotNull final HealthCheck check) {
        super(checkType);
        this.check = check;
    }

    @NotNull
    public HealthCheck getCheck() {
        return check;
    }
}
