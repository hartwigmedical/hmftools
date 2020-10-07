package com.hartwig.hmftools.serve.actionability;

import static org.junit.Assert.assertNotNull;

import org.junit.Test;

public class ActionableEventFactoryTest {

    @Test
    public void canResolveFromFileTypeForValid() {
        assertNotNull(ActionableEventFactory.directionFromFileValue("Responsive"));
    }

    @Test (expected = IllegalStateException.class)
    public void exceptionOnInvalidDirectionFileValue() {
        ActionableEventFactory.directionFromFileValue("Whatever");
    }
}