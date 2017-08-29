package com.hartwig.hmftools.healthchecker;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.healthchecker.runners.CheckType;
import com.hartwig.hmftools.healthchecker.runners.HealthChecker;
import com.hartwig.hmftools.healthchecker.runners.SomaticVariantsChecker;

import org.junit.Test;

public class HealthChecksFlyweightTest {

    @Test
    public void canGenerateExistingChecker() throws NotFoundException {
        final HealthChecksFlyweight healthChecksFlyweight = HealthChecksFlyweight.getInstance();
        final HealthChecker somaticChecker = healthChecksFlyweight.getChecker(CheckType.SOMATIC_VARIANTS.toString());
        assertTrue(somaticChecker instanceof SomaticVariantsChecker);
    }

    @Test
    public void canCreateAllCheckers() throws NotFoundException {
        final HealthChecksFlyweight healthChecksFlyweight = HealthChecksFlyweight.getInstance();
        for (final CheckType type : CheckType.values()) {
            assertNotNull(healthChecksFlyweight.getChecker(type.toString()));
        }
    }

    @Test(expected = NotFoundException.class)
    public void yieldExceptionOnUnknownChecker() throws NotFoundException {
        final HealthChecksFlyweight healthChecksFlyweight = HealthChecksFlyweight.getInstance();
        healthChecksFlyweight.getChecker("bla");
    }
}
