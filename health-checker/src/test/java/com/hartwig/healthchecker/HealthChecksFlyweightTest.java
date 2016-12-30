package com.hartwig.healthchecker;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.healthchecker.runners.MappingChecker;
import com.hartwig.healthchecker.runners.CheckType;
import com.hartwig.healthchecker.runners.HealthChecker;
import com.hartwig.healthchecker.exception.NotFoundException;

import org.junit.Test;

public class HealthChecksFlyweightTest {

    @Test
    public void canGenerateExistingChecker() throws NotFoundException {
        final HealthChecksFlyweight healthChecksFlyweight = HealthChecksFlyweight.getInstance();
        final HealthChecker mappingChecker = healthChecksFlyweight.getChecker(CheckType.MAPPING.toString());
        assertTrue(mappingChecker instanceof MappingChecker);
    }

    @Test
    public void canCreateAllCheckers() throws NotFoundException {
        final HealthChecksFlyweight healthChecksFlyweight = HealthChecksFlyweight.getInstance();
        for (CheckType type : CheckType.values()) {
            assertNotNull(healthChecksFlyweight.getChecker(type.toString()));
        }
    }

    @Test(expected = NotFoundException.class)
    public void yieldExceptionOnUnknownChecker() throws NotFoundException {
        final HealthChecksFlyweight healthChecksFlyweight = HealthChecksFlyweight.getInstance();
        healthChecksFlyweight.getChecker("bla");
    }
}
