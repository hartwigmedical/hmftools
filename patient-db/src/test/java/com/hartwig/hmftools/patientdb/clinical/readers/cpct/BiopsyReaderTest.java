package com.hartwig.hmftools.patientdb.clinical.readers.cpct;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.patientdb.clinical.curators.CuratorTestFactory;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.CuratedPrimaryTumor;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableCuratedPrimaryTumor;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class BiopsyReaderTest {

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @Test
    public void canFilterOutEmptyAndDuplicateForms() {
        BiopsyReader biopsyReader = new BiopsyReader(CuratorTestFactory.biopsySiteCurator());
        List<BiopsyData> biopsies = biopsyReader.read(buildTestPatient(), noCuration());
        assertEquals(1, biopsies.size());
        assertNotNull(biopsies.get(0).date());
        assertNotNull(biopsies.get(0).biopsyTaken());
    }

    @Test
    public void canReadCpctPatientBiopsies() {
        BiopsyReader biopsyReader = new BiopsyReader(CuratorTestFactory.biopsySiteCurator());
        List<BiopsyData> biopsies = biopsyReader.read(buildTestPatient(), noCuration());

        assertEquals(1, biopsies.size());
        assertEquals("body", biopsies.get(0).site());
        assertEquals("liver", biopsies.get(0).location());
        assertEquals(LocalDate.parse("2017-11-10", DATE_FORMATTER), biopsies.get(0).date());
    }

    @NotNull
    private static EcrfPatient buildTestPatient() {
        String patient = "dummy";

        EcrfItemGroup biopsy = new EcrfItemGroup();
        biopsy.addItem(BiopsyReader.FIELD_BIOPSY_TAKEN, "Yes");
        biopsy.addItem(BiopsyReader.FIELD_BIOPSY_EVALUABLE, "Yes");

        // Create initial biopsy.
        EcrfItemGroup biopsies1 = new EcrfItemGroup();
        biopsies1.addItem(BiopsyReader.FIELD_BIOPSY_DATE, "2017-11-10");
        biopsies1.addItem(BiopsyReader.FIELD_SITE, "other");
        biopsies1.addItem(BiopsyReader.FIELD_SITE_OTHER, "body");
        biopsies1.addItem(BiopsyReader.FIELD_LOCATION, "liver");

        EcrfForm form1 = new EcrfForm(FormStatus.undefined());
        form1.addItemGroup(BiopsyReader.ITEMGROUP_BIOPSY, biopsy);
        form1.addItemGroup(BiopsyReader.ITEMGROUP_BIOPSIES, biopsies1);

        // Create empty 2nd biopsy
        EcrfItemGroup biopsies2 = new EcrfItemGroup();
        EcrfForm form2 = new EcrfForm(FormStatus.undefined());
        form2.addItemGroup(BiopsyReader.ITEMGROUP_BIOPSY, biopsy);
        form2.addItemGroup(BiopsyReader.ITEMGROUP_BIOPSIES, biopsies2);

        // Create duplicate 3rd biopsy
        EcrfItemGroup biopsies3 = new EcrfItemGroup();
        biopsies3.addItem(BiopsyReader.FIELD_BIOPSY_DATE, "2017-11-10");
        biopsies3.addItem(BiopsyReader.FIELD_SITE, "other");
        biopsies3.addItem(BiopsyReader.FIELD_SITE_OTHER, "body");
        biopsies3.addItem(BiopsyReader.FIELD_LOCATION, "liver");

        EcrfForm form3 = new EcrfForm(FormStatus.undefined());
        form3.addItemGroup(BiopsyReader.ITEMGROUP_BIOPSY, biopsy);
        form3.addItemGroup(BiopsyReader.ITEMGROUP_BIOPSIES, biopsies3);

        EcrfStudyEvent studyEvent = new EcrfStudyEvent();
        studyEvent.addForm(BiopsyReader.FORM_BIOPS, form1);
        studyEvent.addForm(BiopsyReader.FORM_BIOPS, form2);
        studyEvent.addForm(BiopsyReader.FORM_BIOPS, form3);

        Map<String, List<EcrfStudyEvent>> studyEvents = Maps.newHashMap();
        studyEvents.put(BiopsyReader.STUDY_BIOPSY, Lists.newArrayList(studyEvent));
        return new EcrfPatient(patient, studyEvents, Lists.newArrayList());
    }

    @NotNull
    private static CuratedPrimaryTumor noCuration() {
        return ImmutableCuratedPrimaryTumor.builder().searchTerm(Strings.EMPTY).isOverridden(false).build();
    }
}