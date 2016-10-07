package com.hartwig.hmftools.ecrfanalyser.reader;

import static org.junit.Assert.assertEquals;

import java.io.FileInputStream;
import java.io.FileNotFoundException;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import com.google.common.io.Resources;

import org.junit.Test;

public class XMLEcrfDatamodelReaderTest {

    private static final String TEST_ECRF = Resources.getResource("tests/datamodel.xml").getPath();

    @Test
    public void canExtractDatamodelFromEcrf() throws FileNotFoundException, XMLStreamException {
        XMLInputFactory factory = XMLInputFactory.newInstance();
        XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(TEST_ECRF));
        XMLEcrfDatamodel datamodel = XMLEcrfDatamodelReader.readXMLDatamodel(reader);

        assertEquals(1, datamodel.studyEvents().size());
        assertEquals(1, datamodel.forms().size());
        assertEquals(1, datamodel.itemGroups().size());
        assertEquals(2, datamodel.items().size());
        assertEquals(2, datamodel.codeLists().size());
    }
}