package com.hartwig.hmftools.common.clinical;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.junit.Test;

public class PatientPrimaryTumorFileTest {

    private static final String BASE_RESOURCE_DIR = Resources.getResource("clinical").getPath();
    private static final String TEST_TSV = BASE_RESOURCE_DIR + File.separator + "patient_primary_tumor.tsv";

    @Test
    public void canConvertStringLists() {
        String string1 = "str1";
        String string2 = "str2";

        String fullString = PatientPrimaryTumorFile.fromStringList(Lists.newArrayList(string1, string2));
        List<String> convertedStringList = PatientPrimaryTumorFile.toStringList(fullString);

        assertEquals(2, convertedStringList.size());
        assertEquals(string1, convertedStringList.get(0));
        assertEquals(string2, convertedStringList.get(1));
    }

    @Test
    public void canReadFileAndConvertLines() throws IOException {
        List<PatientPrimaryTumor> patientPrimaryTumors = PatientPrimaryTumorFile.read(TEST_TSV);

        assertEquals(3, patientPrimaryTumors.size());

        List<PatientPrimaryTumor> convertedPrimaryTumors =
                PatientPrimaryTumorFile.fromLines(PatientPrimaryTumorFile.toLines(patientPrimaryTumors));

        assertEquals(3, convertedPrimaryTumors.size());
        for (int i = 0; i < convertedPrimaryTumors.size(); i++) {
            assertEquals(patientPrimaryTumors.get(i), convertedPrimaryTumors.get(i));
        }
    }
}