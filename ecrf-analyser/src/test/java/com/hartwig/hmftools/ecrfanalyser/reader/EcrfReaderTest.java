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

import org.junit.Test;

public class EcrfReaderTest {

    private static final String TEST_ECRF = Resources.getResource("tests/datamodel.xml").getPath();

    @Test
    public void canExtractODMFromEcrf() throws FileNotFoundException, XMLStreamException {
        XMLInputFactory factory = XMLInputFactory.newInstance();
        XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(TEST_ECRF));

        ODMContainer container = EcrfReader.extractODM(reader);
        assertEquals(2, container.itemDefs().size());
        assertEquals(2, container.codeLists().size());
    }

    @Test
    public void canConvertODMContainerToEcrfFields() {
        String category = "blaCategory";
        String fieldName = "blaName";
        String description = "bla";

        String codeListOID = "list";
        String option1 = "x";
        String option2 = "y";

        List<ItemDef> itemDefs = Lists.newArrayList(
                new ItemDef(ItemDefTestFunctions.toOID(category, fieldName), description, codeListOID));
        Map<Integer, String> codeListItems = Maps.newHashMap();
        codeListItems.put(1, option1);
        codeListItems.put(2, option2);

        List<CodeList> codeLists = Lists.newArrayList(new CodeList(codeListOID, codeListItems));
        ODMContainer matchingContainer = new ODMContainer(itemDefs, codeLists);

        List<EcrfField> fields = EcrfReader.odmToEcrfFields(matchingContainer);

        assertEquals(1, fields.size());
        EcrfField field = fields.get(0);
        assertEquals(category, field.category());
        assertEquals(fieldName, field.fieldName());
        assertEquals(description, field.description());

        Map<Integer, String> values = field.values();
        assertEquals(option1, values.get(1));
        assertEquals(option2, values.get(2));
    }
}