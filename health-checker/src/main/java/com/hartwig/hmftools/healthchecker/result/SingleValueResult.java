package com.hartwig.hmftools.healthchecker.result;

import com.hartwig.hmftools.healthchecker.runners.CheckType;
import com.hartwig.hmftools.healthchecker.runners.HealthCheck;

import org.jetbrains.annotations.NotNull;

public class SingleValueResult extends AbstractResult {

    @NotNull
    private final HealthCheck check;

    public SingleValueResult(@NotNull final CheckType checkType, @NotNull final HealthCheck check) {
        super(checkType);
        this.check = check;
    }

    @NotNull
    public HealthCheck check() {
        return check;
    }
}
