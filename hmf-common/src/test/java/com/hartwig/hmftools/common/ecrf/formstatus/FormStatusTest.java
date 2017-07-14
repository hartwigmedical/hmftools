package com.hartwig.hmftools.common.ecrf.formstatus;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.junit.Test;

public class FormStatusTest {

    private static final String BASE_RESOURCE_DIR = Resources.getResource("ecrf").getPath();
    private static final String TEST_FILE = BASE_RESOURCE_DIR + File.separator + "formstatus" + File.separator + "formstatus.csv";

    @Test
    public void canLoadFromCsv() throws IOException, EmptyFileException {
        final FormStatusModel formStatusModel = FormStatus.buildModelFromCsv(TEST_FILE);

        final Map<FormStatusKey, FormStatusData> formStatuses = formStatusModel.formStatuses();
        assertEquals(3, formStatuses.size());
        final FormStatusKey key1 = new ImmutableFormStatusKey("CPCT02000001", "Anti Coagulants (0)", "0", "Anti Coagulants (0)", "0");
        final FormStatusKey key2 =
                new ImmutableFormStatusKey("CPCT02000002", "Death Page (0)", "0", "Neoadjuvant treatment, recurrence and survival (0)",
                        "0");
        final FormStatusKey key3 = new ImmutableFormStatusKey("CPCT02000004", "Eligibility Screening (0)", "0", "BASELINE (0)", "0");

        assertTrue(formStatuses.containsKey(key1));
        assertTrue(formStatuses.containsKey(key2));
        assertTrue(formStatuses.containsKey(key3));

        assertEquals("TRUE", formStatuses.get(key1).locked());
        assertEquals("TRUE", formStatuses.get(key2).locked());
        assertEquals("FALSE", formStatuses.get(key3).locked());

        assertEquals("1", formStatuses.get(key1).dataStatus());
        assertEquals("3", formStatuses.get(key2).dataStatus());
        assertEquals("4", formStatuses.get(key3).dataStatus());
    }
}
