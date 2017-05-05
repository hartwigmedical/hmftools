package com.hartwig.hmftools.common.ecrf.reader;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfField;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class XMLPatientReaderTest {

    private static final String BASE_RESOURCE_DIR = Resources.getResource("ecrf").getPath();
    private static final String PATIENTS_TEST =
            BASE_RESOURCE_DIR + File.separator + "tests" + File.separator + "patients.xml";

    private static final String PATIENT_1 = "CPCT02020202";
    private static final String PATIENT_2 = "CPCT03030303";

    @Test
    public void canReadPatients() throws FileNotFoundException, XMLStreamException {
        final XMLInputFactory factory = XMLInputFactory.newInstance();
        final XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(PATIENTS_TEST));

        final String study = "SE.Study";
        final String form = "FRM.Form";
        final String itemGroup = "GRP.ItemGroup";
        final Map<Integer, String> field1Values = Maps.newHashMap();
        field1Values.put(1, "one");
        field1Values.put(2, "two");
        field1Values.put(3, "three");
        final EcrfField field1 = new EcrfField(study, form, itemGroup, "GRP.ItemGroup.field1", "", field1Values);
        final EcrfField field2 = new EcrfField(study, form, itemGroup, "GRP.ItemGroup.field2", "", Maps.newHashMap());
        final EcrfField birthDate = new EcrfField(study, form, itemGroup, "GRP.ItemGroup.BIRTHDTC", "",
                Maps.newHashMap());

        final List<EcrfPatient> patients = XMLPatientReader.readPatients(reader,
                Lists.newArrayList(field1, field2, birthDate));

        assertEquals(3, patients.size());

        assertEquals(3, patients.get(0).fields().size());
        assertEquals(PATIENT_1, patients.get(0).patientId());
        verifyFirstFieldValue("one", patients.get(0).fieldValuesByEcrfField(field1));
        verifyFirstFieldValue("hi", patients.get(0).fieldValuesByEcrfField(field2));
        verifyFirstFieldValue("2016-01-01", patients.get(0).fieldValuesByEcrfField(birthDate));
        assertEquals(1, patients.get(0).studyEventsPerOID().size());
        assertEquals(1, patients.get(0).studyEventsPerOID().get(study).size());
        assertEquals(1, patients.get(0).studyEventsPerOID().get(study).get(0).formsPerOID().get(form).size());
        assertEquals(1, patients.get(0).studyEventsPerOID().get(study).get(0).formsPerOID().get(form).get(
                0).itemGroupsPerOID().get(itemGroup).size());
        assertEquals(3, patients.get(0).studyEventsPerOID().get(study).get(0).formsPerOID().get(form).get(
                0).itemGroupsPerOID().get(itemGroup).get(0).itemsPerOID().size());

        assertEquals(3, patients.get(1).fields().size());
        assertEquals(PATIENT_2, patients.get(1).patientId());
        verifyFirstFieldValue("two", patients.get(1).fieldValuesByEcrfField(field1));
        verifyFirstFieldValue("hi there", patients.get(1).fieldValuesByEcrfField(field2));
        verifyFirstFieldValue("2016-01-01", patients.get(1).fieldValuesByEcrfField(birthDate));
        assertEquals(1, patients.get(1).studyEventsPerOID().size());
        assertEquals(1, patients.get(1).studyEventsPerOID().get(study).size());
        assertEquals(1, patients.get(1).studyEventsPerOID().get(study).get(0).formsPerOID().get(form).size());
        assertEquals(1, patients.get(1).studyEventsPerOID().get(study).get(0).formsPerOID().get(form).get(
                0).itemGroupsPerOID().get(itemGroup).size());
        assertEquals(3, patients.get(1).studyEventsPerOID().get(study).get(0).formsPerOID().get(form).get(
                0).itemGroupsPerOID().get(itemGroup).get(0).itemsPerOID().size());
    }

    @Test
    public void determinesEmptyItemGroupAndForm() throws FileNotFoundException, XMLStreamException {
        final XMLInputFactory factory = XMLInputFactory.newInstance();
        final XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(PATIENTS_TEST));

        final String study = "SE.Study";
        final String form = "FRM.Form";
        final String itemGroup = "GRP.ItemGroup";
        final Map<Integer, String> field1Values = Maps.newHashMap();
        field1Values.put(1, "one");
        field1Values.put(2, "two");
        field1Values.put(3, "three");
        final EcrfField field1 = new EcrfField(study, form, itemGroup, "GRP.ItemGroup.field1", "", field1Values);
        final EcrfField field2 = new EcrfField(study, form, itemGroup, "GRP.ItemGroup.field2", "", Maps.newHashMap());
        final EcrfField birthDate = new EcrfField(study, form, itemGroup, "GRP.ItemGroup.BIRTHDTC", "",
                Maps.newHashMap());

        final List<EcrfPatient> patients = XMLPatientReader.readPatients(reader,
                Lists.newArrayList(field1, field2, birthDate));

        assertEquals(3, patients.size());

        assertEquals(1, patients.get(2).studyEventsPerOID().size());
        assertEquals(1, patients.get(2).studyEventsPerOID().get(study).size());
        assertEquals(2, patients.get(2).studyEventsPerOID().get(study).get(0).formsPerOID().get(form).size());
        assertEquals(1, patients.get(2).studyEventsPerOID().get(study).get(0).formsPerOID().get(form).get(
                0).itemGroupsPerOID().get(itemGroup).size());
        assertEquals(true, patients.get(2).studyEventsPerOID().get(study).get(0).formsPerOID().get(form).get(
                0).itemGroupsPerOID().get(itemGroup).get(0).isEmpty());
        assertEquals(true,
                patients.get(2).studyEventsPerOID().get(study).get(0).formsPerOID().get(form).get(0).isEmpty());
        assertEquals(2, patients.get(2).studyEventsPerOID().get(study).get(0).formsPerOID().get(form).get(
                1).itemGroupsPerOID().get(itemGroup).size());
        assertEquals(true, patients.get(2).studyEventsPerOID().get(study).get(0).formsPerOID().get(form).get(
                1).itemGroupsPerOID().get(itemGroup).get(0).isEmpty());
        System.out.println(patients.get(2).studyEventsPerOID().get(study).get(0).formsPerOID().get(form).get(
                1).itemGroupsPerOID().get(itemGroup).get(1));
        assertEquals(false, patients.get(2).studyEventsPerOID().get(study).get(0).formsPerOID().get(form).get(
                1).itemGroupsPerOID().get(itemGroup).get(1).isEmpty());
        assertEquals(false,
                patients.get(2).studyEventsPerOID().get(study).get(0).formsPerOID().get(form).get(1).isEmpty());
    }

    private static void verifyFirstFieldValue(@NotNull final String expected, @Nullable final List<String> values) {
        assert values != null && values.size() == 1;
        assertEquals(expected, values.get(0));
    }
}