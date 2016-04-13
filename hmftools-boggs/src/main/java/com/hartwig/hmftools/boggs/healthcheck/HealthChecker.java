package com.hartwig.hmftools.boggs.healthcheck;

import com.hartwig.hmftools.boggs.PatientData;
import org.jetbrains.annotations.NotNull;

public interface HealthChecker {

    boolean isHealthy(@NotNull PatientData patient);
}
