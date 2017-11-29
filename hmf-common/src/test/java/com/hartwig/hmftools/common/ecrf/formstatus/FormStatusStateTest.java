package com.hartwig.hmftools.common.ecrf.formstatus;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class FormStatusStateTest {

    @Test
    public void canDetermineBest() {
        assertEquals(FormStatusState.VERIFIED,
                FormStatusState.best(FormStatusState.UNKNOWN, FormStatusState.SUBMITTED_WITH_MISSING, FormStatusState.VERIFIED,
                        FormStatusState.SAVED));

        assertEquals(FormStatusState.SUBMITTED, FormStatusState.best(FormStatusState.SUBMITTED));
    }
}