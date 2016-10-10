package com.hartwig.hmftools.ecrfanalyser.reader;

import static org.junit.Assert.assertEquals;

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
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfField;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfPatient;

import org.junit.Test;

public class XMLPatientReaderTest {

    private static final String TEST_ECRF = Resources.getResource("tests/patients.xml").getPath();

    private static final String PATIENT_1 = "CPCT02020202";
    private static final String PATIENT_2 = "CPCT03030303";

    @Test
    public void canReadPatients() throws FileNotFoundException, XMLStreamException {
        XMLInputFactory factory = XMLInputFactory.newInstance();
        XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(TEST_ECRF));

        String study = "SE.Study";
        String form = "FRM.Form";
        String itemGroup = "GRP.ItemGroup";
        Map<Integer, String> field1Values = Maps.newHashMap();
        field1Values.put(1, "one");
        field1Values.put(2, "two");
        field1Values.put(3, "three");
        EcrfField field1 = new EcrfField(study, form, itemGroup, "GRP.ItemGroup.field1", "", field1Values);
        EcrfField field2 = new EcrfField(study, form, itemGroup, "GRP.ItemGroup.field2", "",
                Maps.<Integer, String>newHashMap());
        EcrfField birthDate = new EcrfField(study, form, itemGroup, "GRP.ItemGroup.BIRTHDTC", "",
                Maps.<Integer, String>newHashMap());

        List<EcrfPatient> patients = XMLPatientReader.readPatients(reader,
                Lists.newArrayList(field1, field2, birthDate));

        assertEquals(2, patients.size());

        assertEquals(3, patients.get(0).fields().size());
        assertEquals(PATIENT_1, patients.get(0).patientId());
        assertEquals("one", patients.get(0).fieldValue(field1));
        assertEquals("hi", patients.get(0).fieldValue(field2));
        assertEquals("2016-01-01", patients.get(0).fieldValue(birthDate));

        assertEquals(3, patients.get(1).fields().size());
        assertEquals(PATIENT_2, patients.get(1).patientId());
        assertEquals("two", patients.get(1).fieldValue(field1));
        assertEquals("hi there", patients.get(1).fieldValue(field2));
        assertEquals("2016-01-01", patients.get(1).fieldValue(birthDate));
    }
}