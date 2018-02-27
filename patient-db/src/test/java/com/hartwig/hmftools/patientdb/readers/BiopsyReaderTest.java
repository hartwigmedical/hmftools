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
import com.hartwig.hmftools.patientdb.data.BiopsyData;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class BiopsyReaderTest {

    @Test
    public void canFilterOutEmptyAndDuplicateForms() {
        List<BiopsyData> biopsies = BiopsyReader.read(buildTestPatient());
        assertEquals(1, biopsies.size());
        assertNotNull(biopsies.get(0).date());
        assertNotNull(biopsies.get(0).biopsyTaken());
    }

    @NotNull
    private static EcrfPatient buildTestPatient() {
        final String patient = "dummy";

        EcrfItemGroup biopsy = new EcrfItemGroup(patient);
        biopsy.addItem(BiopsyReader.FIELD_BIOPSY_TAKEN, "Yes");
        biopsy.addItem(BiopsyReader.FIELD_BIOPSY_EVALUABLE, "Yes");

        // KODU: Create initial biopsy.
        EcrfItemGroup biopsies1 = new EcrfItemGroup(patient);
        biopsies1.addItem(BiopsyReader.FIELD_BIOPSY_DATE, "2017-11-10");
        biopsies1.addItem(BiopsyReader.FIELD_LOCATION, "liver");
        biopsies1.addItem(BiopsyReader.FIELD_SITE, "other");
        biopsies1.addItem(BiopsyReader.FIELD_SITE_OTHER, "body");

        EcrfForm form1 = new EcrfForm(patient, FormStatusState.SAVED, true);
        form1.addItemGroup(BiopsyReader.ITEMGROUP_BIOPSY, biopsy);
        form1.addItemGroup(BiopsyReader.ITEMGROUP_BIOPSIES, biopsies1);

        // KODU: Create empty 2nd biopsy
        EcrfItemGroup biopsies2 = new EcrfItemGroup(patient);
        EcrfForm form2 = new EcrfForm(patient, FormStatusState.SAVED, true);
        form2.addItemGroup(BiopsyReader.ITEMGROUP_BIOPSY, biopsy);
        form2.addItemGroup(BiopsyReader.ITEMGROUP_BIOPSIES, biopsies2);

        // KODU: Create duplicate 3rd biopsy
        EcrfItemGroup biopsies3 = new EcrfItemGroup(patient);
        biopsies3.addItem(BiopsyReader.FIELD_BIOPSY_DATE, "2017-11-10");
        biopsies3.addItem(BiopsyReader.FIELD_LOCATION, "liver");
        biopsies3.addItem(BiopsyReader.FIELD_SITE, "other");
        biopsies3.addItem(BiopsyReader.FIELD_SITE_OTHER, "body");

        EcrfForm form3 = new EcrfForm(patient, FormStatusState.SAVED, true);
        form3.addItemGroup(BiopsyReader.ITEMGROUP_BIOPSY, biopsy);
        form3.addItemGroup(BiopsyReader.ITEMGROUP_BIOPSIES, biopsies3);

        EcrfStudyEvent studyEvent = new EcrfStudyEvent(patient);
        studyEvent.addForm(BiopsyReader.FORM_BIOPS, form1);
        studyEvent.addForm(BiopsyReader.FORM_BIOPS, form2);
        studyEvent.addForm(BiopsyReader.FORM_BIOPS, form3);

        Map<String, List<EcrfStudyEvent>> studyEvents = Maps.newHashMap();
        studyEvents.put(BiopsyReader.STUDY_BIOPSY, Lists.newArrayList(studyEvent));
        return new EcrfPatient(patient, studyEvents, Lists.newArrayList());
    }
}