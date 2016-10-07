package com.hartwig.hmftools.ecrfanalyser.reader;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfField;

import org.junit.Test;

public class XMLEcrfDatamodelToEcrfFieldsTest {

    @Test
    public void canConvertXMLObjectContainerToEcrfFields() {
        String fieldName = "blaCategory.blaName";
        String description = "bla";

        String codeListOID = "list";
        String option1 = "x";
        String option2 = "y";

        List<Item> items = Lists.newArrayList(new Item(OIDFunctions.toOID(fieldName), description, codeListOID));
        Map<Integer, String> codeListItems = Maps.newHashMap();
        codeListItems.put(1, option1);
        codeListItems.put(2, option2);

        List<CodeList> codeLists = Lists.newArrayList(new CodeList(codeListOID, codeListItems));
        XMLEcrfDatamodel matchingContainer = new XMLEcrfDatamodel(Lists.<StudyEvent>newArrayList(),
                Lists.<Form>newArrayList(), Lists.<ItemGroup>newArrayList(), items, codeLists);

        List<EcrfField> fields = XMLEcrfDatamodelToEcrfFields.convert(matchingContainer);

        assertEquals(1, fields.size());
        EcrfField field = fields.get(0);
        assertEquals(fieldName, field.name());
        assertEquals(description, field.description());

        Map<Integer, String> values = field.values();
        assertEquals(option1, values.get(1));
        assertEquals(option2, values.get(2));
    }
}