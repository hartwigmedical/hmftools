package com.hartwig.hmftools.boggs.healthcheck;

import com.hartwig.hmftools.boggs.PatientData;
import org.jetbrains.annotations.NotNull;

public class MappingHealthChecker implements HealthChecker {

    public boolean isHealthy(@NotNull PatientData patient) {
        return true;
    }
}
