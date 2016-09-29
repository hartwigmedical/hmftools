package com.hartwig.hmftools.ecrfanalyser.reader;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.List;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import com.google.common.io.Resources;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfField;

import org.junit.Test;

public class EcrfReaderTest {


    // KODU: ItemDefs come in 2 flavours: with or without CodeListRef, depending on the DataType.
//    <ItemDef OID="whatever" Name="whatever" DataType="text" Length="whatever" SASFieldName="whatever" />
//    <ItemDef OID="whatever" Name="whatever" DataType="integer" Length="whatever" SASFieldName="whatever">
//        <CodeListRef CodeListOID="whatever" />
//    </ItemDef>

    private static final String TEST_ECRF = Resources.getResource("tests/datamodel.xml").getPath();

    @Test
    public void canExtractODMFromEcrf() throws FileNotFoundException, XMLStreamException {
        XMLInputFactory factory = XMLInputFactory.newInstance();
        XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(TEST_ECRF));

        ODMContainer container = EcrfReader.extractODM(reader);
//        for (ItemDef itemDef : container.itemDefs()) {
//            System.out.pri
//        }
        System.out.println(container.itemDefs());
        System.out.println(container.codeLists());
    }
}