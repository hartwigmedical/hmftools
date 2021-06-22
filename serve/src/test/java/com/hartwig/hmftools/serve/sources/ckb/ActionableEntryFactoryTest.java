package com.hartwig.hmftools.serve.sources.ckb;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import org.junit.Test;

public class ActionableEntryFactoryTest {

    @Test
    public void canExtractAndMapDoid() {
        assertNull(ActionableEntryFactory.extractDoid(null));
        assertNull(ActionableEntryFactory.extractDoid("not a doid"));

        assertEquals("0060463", ActionableEntryFactory.extractDoid("DOID:0060463"));
        assertEquals("162", ActionableEntryFactory.extractDoid("JAX:10000003"));
        assertNull(ActionableEntryFactory.extractDoid("JAX:10000004"));
    }
}
