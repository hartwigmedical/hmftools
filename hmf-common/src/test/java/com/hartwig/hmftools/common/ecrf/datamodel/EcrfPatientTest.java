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
        final String patientId = "patient";
        final EcrfField field = new ImmutableEcrfField("1", "2", "3", "4", "5", Maps.newHashMap());
        final EcrfItemGroup itemGroup = new EcrfItemGroup(patientId);
        itemGroup.itemsPerOID().put("4", Lists.newArrayList("value1", "value2"));
        final EcrfForm form = new EcrfForm(patientId, "", "");
        form.itemGroupsPerOID().put("3", Lists.newArrayList(itemGroup));
        final EcrfStudyEvent study = new EcrfStudyEvent(patientId);
        study.formsPerOID().put("2", Lists.newArrayList(form));
        final Map<String, List<EcrfStudyEvent>> studyEvents = Maps.newHashMap();
        studyEvents.put("1", Lists.newArrayList(study));
        final EcrfPatient patient = new EcrfPatient(patientId, studyEvents);

        assertNull(patient.fieldValuesByName("does not exist"));
        final List<String> values = patient.fieldValuesByName(field.name());
        assertNotNull(values);
        assertEquals(2, values.size());
    }
}