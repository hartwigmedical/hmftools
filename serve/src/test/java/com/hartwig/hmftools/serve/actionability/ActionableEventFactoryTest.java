package com.hartwig.hmftools.serve.actionability;

import static org.junit.Assert.assertNotNull;

import org.junit.Test;

public class ActionableEventFactoryTest {

    @Test
    public void canResolveFromDisplayForValid() {
        assertNotNull(ActionableEventFactory.directionFromFileValue("Responsive"));
        assertNotNull(ActionableEventFactory.sourceFromFileValue("CIViC"));
    }

    @Test (expected = IllegalStateException.class)
    public void exceptionOnInvalidDirectionFileValue() {
        ActionableEventFactory.directionFromFileValue("Whatever");
    }

    @Test (expected = IllegalStateException.class)
    public void exceptionOnInvalidSourceFileValue() {
        ActionableEventFactory.sourceFromFileValue("Whatever");
    }
}