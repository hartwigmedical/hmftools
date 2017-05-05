package com.hartwig.hmftools.common.ecrf.datamodel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.junit.Test;

public class EcrfPatientTest {

    @Test
    public void canFindValuesById() {
        final EcrfField field = new EcrfField("1", "2", "3", "4", "5", Maps.newHashMap());
        final Map<EcrfField, List<String>> fieldValues = Maps.newHashMap();
        fieldValues.put(field, Lists.newArrayList("value1", "value2"));
        final EcrfPatient patient = new EcrfPatient("patient", fieldValues, Maps.newHashMap());

        assertNull(patient.fieldValuesByName("does not exist"));
        final List<String> values = patient.fieldValuesByName(field.name());
        assertNotNull(values);
        assertEquals(2, values.size());
    }

}