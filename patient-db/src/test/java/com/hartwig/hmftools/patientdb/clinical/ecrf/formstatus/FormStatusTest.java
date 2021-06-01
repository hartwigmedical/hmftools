package com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class FormStatusTest {

    @Test
    public void canMerge() {
        FormStatus status1 = ImmutableFormStatus.builder().state(FormStatusState.SAVED).locked(true).build();
        FormStatus status2 = ImmutableFormStatus.builder().state(FormStatusState.SUBMITTED).locked(true).build();
        FormStatus status3 = ImmutableFormStatus.builder().state(FormStatusState.VERIFIED).locked(false).build();

        FormStatus merged1 = FormStatus.merge(status1);
        assertEquals(FormStatusState.SAVED, merged1.state());
        assertTrue(merged1.locked());

        FormStatus merged2 = FormStatus.merge(status1, status2);
        assertEquals(FormStatusState.SUBMITTED, merged2.state());
        assertTrue(merged2.locked());

        FormStatus merged3 = FormStatus.merge(status1, status2, status3);
        assertEquals(FormStatusState.VERIFIED, merged3.state());
        assertFalse(merged3.locked());
    }
}