package com.hartwig.hmftools.patientdb.clinical.ecrf.reader;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import com.google.common.io.Resources;

import org.junit.Test;

public class XMLEcrfDatamodelReaderTest {

    private static final String BASE_RESOURCE_DIR = Resources.getResource("ecrf").getPath();
    private static final String DATAMODEL_TEST = BASE_RESOURCE_DIR + File.separator + "tests" + File.separator + "datamodel.xml";

    @Test
    public void canExtractDatamodelFromEcrf() throws FileNotFoundException, XMLStreamException {
        XMLInputFactory factory = XMLInputFactory.newInstance();
        XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(DATAMODEL_TEST));
        XMLEcrfDatamodel datamodel = XMLEcrfDatamodelReader.readXMLDatamodel(reader);

        assertEquals(1, datamodel.studyEvents().size());
        assertEquals(1, datamodel.forms().size());
        assertEquals(1, datamodel.itemGroups().size());
        assertEquals(2, datamodel.items().size());
        assertEquals(2, datamodel.codeLists().size());
    }
}