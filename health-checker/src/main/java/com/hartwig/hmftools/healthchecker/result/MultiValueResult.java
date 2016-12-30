package com.hartwig.hmftools.healthchecker.result;

import java.util.List;

import com.hartwig.hmftools.healthchecker.runners.CheckType;
import com.hartwig.hmftools.healthchecker.runners.checks.HealthCheck;

import org.jetbrains.annotations.NotNull;

public class MultiValueResult extends AbstractResult {

    private static final long serialVersionUID = 7396913735152649393L;

    @NotNull
    private final List<HealthCheck> checks;

    public MultiValueResult(@NotNull final CheckType checkType, @NotNull final List<HealthCheck> checks) {
        super(checkType);
        this.checks = checks;
    }

    @NotNull
    public List<HealthCheck> getChecks() {
        return checks;
    }
}
