package com.hartwig.hmftools.patientdb.readers;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class BiopsyTreatmentResponseReaderTest {

    @Test
    public void canFilterOutEmptyForms() {
        List<BiopsyTreatmentResponseData> responses = BiopsyTreatmentResponseReader.read(buildTestPatient());
        assertEquals(1, responses.size());
        assertNotNull(responses.get(0).measurementDone());
        assertNotNull(responses.get(0).response());
    }

    @NotNull
    private static EcrfPatient buildTestPatient() {
        final String patient = "dummy";

        EcrfItemGroup response1 = new EcrfItemGroup(patient);
        response1.addItem(BiopsyTreatmentResponseReader.FIELD_MEASUREMENT_DONE, "Yes");
        response1.addItem(BiopsyTreatmentResponseReader.FIELD_RESPONSE, "PR");

        EcrfForm form1 = new EcrfForm(patient, FormStatusState.SAVED, true);
        form1.addItemGroup(BiopsyTreatmentResponseReader.ITEMGROUP_TUMOR_MEASUREMENT, response1);

        EcrfItemGroup response2 = new EcrfItemGroup(patient);

        EcrfForm form2 = new EcrfForm(patient, FormStatusState.SAVED, true);
        form1.addItemGroup(BiopsyTreatmentResponseReader.ITEMGROUP_TUMOR_MEASUREMENT, response2);

        EcrfStudyEvent studyEvent = new EcrfStudyEvent(patient);
        studyEvent.addForm(BiopsyTreatmentResponseReader.FORM_TUMOR_MEASUREMENT, form1);
        studyEvent.addForm(BiopsyTreatmentResponseReader.FORM_TUMOR_MEASUREMENT, form2);

        Map<String, List<EcrfStudyEvent>> studyEvents = Maps.newHashMap();
        studyEvents.put(BiopsyTreatmentResponseReader.STUDY_TREATMENT, Lists.newArrayList(studyEvent));
        return new EcrfPatient(patient, studyEvents, Lists.newArrayList());
    }
}