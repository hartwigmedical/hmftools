package com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.Map;

import com.google.common.io.Resources;

import org.junit.Test;

public class FormStatusReaderTest {

    private static final String TEST_FILE = Resources.getResource("ecrf/test_cpct_formstatus.csv").getPath();

    @Test
    public void canLoadFromCsv() throws IOException {
        FormStatusModel formStatusModel = FormStatusReader.buildModelFromCsv(TEST_FILE);

        Map<FormStatusKey, FormStatus> formStatuses = formStatusModel.formStatuses();
        assertEquals(3, formStatuses.size());

        FormStatusKey key1 = new ImmutableFormStatusKey("CPCT02000001", "Anti Coagulants", "0", "Anti Coagulants", "0");
        FormStatusKey key2 =
                new ImmutableFormStatusKey("CPCT02000002", "Death Page", "0", "Neoadjuvant treatment, recurrence and survival", "0");
        FormStatusKey key3 = new ImmutableFormStatusKey("CPCT02000004", "Eligibility Screening", "0", "BASELINE", "0");

        assertTrue(formStatuses.containsKey(key1));
        assertTrue(formStatuses.containsKey(key2));
        assertTrue(formStatuses.containsKey(key3));

        assertTrue(formStatuses.get(key1).locked());
        assertTrue(formStatuses.get(key2).locked());
        assertFalse(formStatuses.get(key3).locked());

        assertEquals(FormStatusState.SUBMITTED, formStatuses.get(key1).state());
        assertEquals(FormStatusState.SUBMITTED_WITH_MISSING, formStatuses.get(key2).state());
        assertEquals(FormStatusState.VERIFIED, formStatuses.get(key3).state());
    }
}
