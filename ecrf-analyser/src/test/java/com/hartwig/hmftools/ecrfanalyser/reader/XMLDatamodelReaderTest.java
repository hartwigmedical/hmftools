package com.hartwig.hmftools.ecrfanalyser.reader;

import static org.junit.Assert.assertEquals;

import java.io.FileInputStream;
import java.io.FileNotFoundException;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import com.google.common.io.Resources;

import org.junit.Test;

public class XMLDatamodelReaderTest {

    private static final String TEST_ECRF = Resources.getResource("tests/datamodel.xml").getPath();

    @Test
    public void canExtractDatamodelFromEcrf() throws FileNotFoundException, XMLStreamException {
        XMLInputFactory factory = XMLInputFactory.newInstance();
        XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(TEST_ECRF));

        XMLEcrfDatamodel container = XMLDatamodelReader.readXMLDatamodel(reader);
        assertEquals(2, container.items().size());
        assertEquals(2, container.codeLists().size());
    }
}