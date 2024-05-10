package com.hartwig.hmftools.peach.effect;

import static org.junit.Assert.assertNull;

import java.util.Collections;

import org.junit.Test;

public class DrugInfoStoreTest
{
    @Test
    public void testMissingInfo()
    {
        DrugInfoStore store = new DrugInfoStore(Collections.emptyList());
        assertNull(store.getPrescriptionInfoUrl("FAKE_GENE", "FAKE_DRUG"));
    }
}
