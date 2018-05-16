package com.hartwig.hmftools.common.ecrf.datamodel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;

import org.junit.Test;

public class EcrfPatientTest {

    @Test
    public void canFindValuesById() {
        final String patientId = "patient";
        final EcrfDatamodelField field = new ImmutableEcrfDatamodelField("1", "2", "3", "4", "5", Maps.newHashMap());
        final EcrfItemGroup itemGroup = new EcrfItemGroup();
        itemGroup.itemsPerOID().put("4", Lists.newArrayList("value1", "value2"));
        final EcrfForm form = new EcrfForm(FormStatus.undefined());
        form.itemGroupsPerOID().put("3", Lists.newArrayList(itemGroup));
        final EcrfStudyEvent study = new EcrfStudyEvent();
        study.formsPerOID().put("2", Lists.newArrayList(form));
        final Map<String, List<EcrfStudyEvent>> studyEvents = Maps.newHashMap();
        studyEvents.put("1", Lists.newArrayList(study));
        final EcrfPatient patient = new EcrfPatient(patientId, studyEvents,
                Lists.newArrayList(EcrfDataField.of(patientId, "1", "1", "2", "1", "3", "1", "4", "value1", "", ""),
                        EcrfDataField.of(patientId, "1", "1", "2", "1", "3", "2", "4", "value1", "", "")));

        assertNull(patient.fieldValuesByName("does not exist"));
        final List<String> valuesByName = patient.fieldValuesByName(field.name());
        assertNotNull(valuesByName);
        assertEquals(2, valuesByName.size());
        final List<String> valuesByField = patient.fieldValuesByEcrfField(field);
        assertNotNull(valuesByField);
        assertEquals(2, valuesByField.size());
    }
}
