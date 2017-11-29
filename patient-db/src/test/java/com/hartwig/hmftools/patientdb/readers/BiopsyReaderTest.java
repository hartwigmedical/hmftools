package com.hartwig.hmftools.patientdb.readers;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.data.BiopsyData;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class BiopsyReaderTest {

    @Test
    public void canFilterOutEmptyForms() {
        List<BiopsyData> biopsies = BiopsyReader.read(buildTestPatient());
        assertEquals(1, biopsies.size());
    }

    @NotNull
    private static EcrfPatient buildTestPatient() {
        final String patient = "dummy";

        EcrfItemGroup biopsy1 = new EcrfItemGroup(patient);
        biopsy1.addItem(BiopsyReader.FIELD_BIOPSY_DATE, "2017-11-10");
        biopsy1.addItem(BiopsyReader.FIELD_LOCATION, "liver");
        biopsy1.addItem(BiopsyReader.FIELD_SITE, "other");
        biopsy1.addItem(BiopsyReader.FIELD_SITE_OTHER, "body");

        EcrfForm form1 = new EcrfForm(patient, Strings.EMPTY, "TRUE");
        form1.addItemGroup(BiopsyReader.ITEMGROUP_BIOPSIES, biopsy1);

        EcrfItemGroup biopsy2 = new EcrfItemGroup(patient);
        EcrfForm form2 = new EcrfForm(patient, Strings.EMPTY, "TRUE");
        form2.addItemGroup(BiopsyReader.ITEMGROUP_BIOPSIES, biopsy2);

        EcrfStudyEvent studyEvent = new EcrfStudyEvent(patient);
        studyEvent.addForm(BiopsyReader.FORM_BIOPS, form1);
        studyEvent.addForm(BiopsyReader.FORM_BIOPS, form2);

        Map<String, List<EcrfStudyEvent>> studyEvents = Maps.newHashMap();
        studyEvents.put(BiopsyReader.STUDY_BIOPSY, Lists.newArrayList(studyEvent));
        return new EcrfPatient(patient, studyEvents, Lists.newArrayList());
    }
}