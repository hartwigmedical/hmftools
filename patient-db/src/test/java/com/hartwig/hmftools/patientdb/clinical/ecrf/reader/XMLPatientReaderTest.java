package com.hartwig.hmftools.patientdb.clinical.ecrf.reader;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

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
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfDataField;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfDatamodelField;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfField;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.ImmutableEcrfDatamodelField;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.ImmutableFormStatusModel;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class XMLPatientReaderTest {

    private static final String PATIENTS_TEST = Resources.getResource("ecrf/patients_example.xml").getPath();

    private static final String PATIENT_1 = "CPCT02020202";
    private static final String PATIENT_2 = "CPCT03030303";

    @Test
    public void canReadPatients() throws FileNotFoundException, XMLStreamException {
        XMLInputFactory factory = XMLInputFactory.newInstance();
        XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(PATIENTS_TEST));

        String studyOID = "SE.Study";
        String formOID = "FRM.Form";
        String itemGroupOID = "GRP.ItemGroup";
        Map<Integer, String> codeListValues = Maps.newHashMap();
        codeListValues.put(1, "one");
        codeListValues.put(2, "two");
        codeListValues.put(3, "three");
        String item1OID = "FLD.ItemGroup.field1";
        String item2OID = "FLD.ItemGroup.field2";
        String birthDateOID = "FLD.ItemGroup.BIRTHDTC";
        String codeListOID = "codeList";
        StudyEvent studyEvent = new ImmutableStudyEvent(studyOID, studyOID, Lists.newArrayList(formOID));
        Form form = new ImmutableForm(formOID, formOID, Lists.newArrayList(itemGroupOID));
        ItemGroup itemGroup = new ImmutableItemGroup(itemGroupOID, itemGroupOID, Lists.newArrayList(item1OID, item2OID, birthDateOID));
        Item item1 = new ImmutableItem(item1OID, item1OID, codeListOID);
        Item item2 = new ImmutableItem(item2OID, item2OID, "");
        Item birthDate = new ImmutableItem(birthDateOID, birthDateOID, "");
        CodeList codeList = new ImmutableCodeList(codeListOID, codeListOID, codeListValues);

        List<EcrfPatient> patients = XMLPatientReader.readPatients(reader,
                XMLEcrfDatamodel.of(Lists.newArrayList(studyEvent),
                        Lists.newArrayList(form),
                        Lists.newArrayList(itemGroup),
                        Lists.newArrayList(item1, item2, birthDate),
                        Lists.newArrayList(codeList)),
                new ImmutableFormStatusModel(Maps.newHashMap()));

        EcrfDatamodelField field1 =
                new ImmutableEcrfDatamodelField(studyOID, formOID, itemGroupOID, "FLD.ItemGroup.field1", "", codeListValues);
        EcrfDatamodelField field2 =
                new ImmutableEcrfDatamodelField(studyOID, formOID, itemGroupOID, "FLD.ItemGroup.field2", "", Maps.newHashMap());
        EcrfDatamodelField birthDateField =
                new ImmutableEcrfDatamodelField(studyOID, formOID, itemGroupOID, "FLD.ItemGroup.BIRTHDTC", "", Maps.newHashMap());

        // @formatter:off
        assertEquals(3, patients.size());
        assertEquals(4, patients.get(0).fields().size());
        assertEquals(PATIENT_1, patients.get(0).patientId());
        verifyFirstFieldValue("one", fieldValues(patients.get(0), field1));
        verifyFirstFieldValue("hi", fieldValues(patients.get(0), field2));
        verifyFirstFieldValue("2016-01-01", fieldValues(patients.get(0), birthDateField));
        verifyFirstFieldValue("one", patients.get(0).studyEventsPerOID(studyOID).get(0).formsPerOID().get(formOID).get(0).itemGroupsPerOID()
                .get(itemGroupOID).get(0).itemsPerOID().get(item1OID));
        verifyFirstFieldValue("hi", patients.get(0).studyEventsPerOID(studyOID).get(0).formsPerOID().get(formOID).get(0).itemGroupsPerOID()
                .get(itemGroupOID).get(0).itemsPerOID().get(item2OID));
        verifyFirstFieldValue("2016-01-01", patients.get(0).studyEventsPerOID(studyOID).get(0).formsPerOID().get(formOID).get(0).itemGroupsPerOID()
                .get(itemGroupOID).get(0).itemsPerOID().get(birthDateOID));
        assertEquals(1, patients.get(0).studyEventsPerOID().size());
        assertEquals(1, patients.get(0).studyEventsPerOID().get(studyOID).size());
        assertEquals(1, patients.get(0).studyEventsPerOID().get(studyOID).get(0).formsPerOID().get(formOID).size());
        assertEquals(1, patients.get(0).studyEventsPerOID().get(studyOID).get(0).formsPerOID().get(formOID).get(
                0).itemGroupsPerOID().get(itemGroupOID).size());
        assertEquals(4, patients.get(0).studyEventsPerOID().get(studyOID).get(0).formsPerOID().get(formOID).get(
                0).itemGroupsPerOID().get(itemGroupOID).get(0).itemsPerOID().size());

        assertEquals(4, patients.get(1).fields().size());
        assertEquals(PATIENT_2, patients.get(1).patientId());
        verifyFirstFieldValue("two", patients.get(1).studyEventsPerOID(studyOID).get(0).formsPerOID().get(formOID).get(0).itemGroupsPerOID()
                .get(itemGroupOID).get(0).itemsPerOID().get(item1OID));
        verifyFirstFieldValue("hi there", patients.get(1).studyEventsPerOID(studyOID).get(0).formsPerOID().get(formOID).get(0).itemGroupsPerOID()
                .get(itemGroupOID).get(0).itemsPerOID().get(item2OID));
        verifyFirstFieldValue("2016-01-01", patients.get(1).studyEventsPerOID(studyOID).get(0).formsPerOID().get(formOID).get(0).itemGroupsPerOID()
                .get(itemGroupOID).get(0).itemsPerOID().get(birthDateOID));
        assertEquals(1, patients.get(1).studyEventsPerOID().size());
        assertEquals(1, patients.get(1).studyEventsPerOID().get(studyOID).size());
        assertEquals(1, patients.get(1).studyEventsPerOID().get(studyOID).get(0).formsPerOID().get(formOID).size());
        assertEquals(1, patients.get(1).studyEventsPerOID().get(studyOID).get(0).formsPerOID().get(formOID).get(
                0).itemGroupsPerOID().get(itemGroupOID).size());
        assertEquals(4, patients.get(1).studyEventsPerOID().get(studyOID).get(0).formsPerOID().get(formOID).get(
                0).itemGroupsPerOID().get(itemGroupOID).get(0).itemsPerOID().size());
        // @formatter:on
    }

    @Test
    public void determinesEmptyItemGroupAndForm() throws FileNotFoundException, XMLStreamException {
        XMLInputFactory factory = XMLInputFactory.newInstance();
        XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(PATIENTS_TEST));

        String studyOID = "SE.Study";
        String formOID = "FRM.Form";
        String itemGroupOID = "GRP.ItemGroup";

        StudyEvent studyEvent = new ImmutableStudyEvent(studyOID, studyOID, Lists.newArrayList(formOID));
        Form form = new ImmutableForm(formOID, formOID, Lists.newArrayList(itemGroupOID));
        List<EcrfPatient> patients = XMLPatientReader.readPatients(reader,
                XMLEcrfDatamodel.of(Lists.newArrayList(studyEvent),
                        Lists.newArrayList(form),
                        Lists.newArrayList(),
                        Lists.newArrayList(),
                        Lists.newArrayList()),
                new ImmutableFormStatusModel(Maps.newHashMap()));

        // @formatter:off
        assertEquals(3, patients.size());
        assertEquals(1, patients.get(2).studyEventsPerOID().size());
        assertEquals(1, patients.get(2).studyEventsPerOID().get(studyOID).size());
        assertEquals(2, patients.get(2).studyEventsPerOID().get(studyOID).get(0).formsPerOID().get(formOID).size());
        assertEquals(1, patients.get(2).studyEventsPerOID().get(studyOID).get(0).formsPerOID().get(formOID).get(
                0).itemGroupsPerOID().get(itemGroupOID).size());
        assertTrue(patients.get(2).studyEventsPerOID().get(studyOID).get(0).formsPerOID().get(formOID).get(
                0).itemGroupsPerOID().get(itemGroupOID).get(0).isEmpty());
        assertTrue(patients.get(2).studyEventsPerOID().get(studyOID).get(0).formsPerOID().get(formOID).get(0).isEmpty());
        assertEquals(2, patients.get(2).studyEventsPerOID().get(studyOID).get(0).formsPerOID().get(formOID).get(
                1).itemGroupsPerOID().get(itemGroupOID).size());
        assertTrue(patients.get(2).studyEventsPerOID().get(studyOID).get(0).formsPerOID().get(formOID).get(
                1).itemGroupsPerOID().get(itemGroupOID).get(0).isEmpty());
        assertFalse(patients.get(2).studyEventsPerOID().get(studyOID).get(0).formsPerOID().get(formOID).get(
                1).itemGroupsPerOID().get(itemGroupOID).get(1).isEmpty());
        assertFalse(patients.get(2).studyEventsPerOID().get(studyOID).get(0).formsPerOID().get(formOID).get(1).isEmpty());
        // @formatter:on
    }

    @NotNull
    private static List<String> fieldValues(@NotNull EcrfPatient patient, @NotNull EcrfField field) {
        List<String> fieldValues = Lists.newArrayList();
        for (EcrfDataField dataField : patient.fields()) {
            if (dataField.name().equals(field.name())) {
                fieldValues.add(dataField.itemValue());
            }
        }
        return fieldValues;
    }

    private static void verifyFirstFieldValue(@NotNull String expected, @Nullable List<String> values) {
        assert values != null && values.size() == 1;
        assertEquals(expected, values.get(0));
    }
}