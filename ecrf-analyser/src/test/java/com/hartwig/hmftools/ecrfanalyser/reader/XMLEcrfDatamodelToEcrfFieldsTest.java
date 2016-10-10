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
        String studyOID = "SE.Study";
        String formOID = "FRM.Form";
        String itemGroupOID = "GRP.ItemGroup";
        String itemOID = "GRP.ItemGroup.Item";
        String description = "bla";

        String codeListOID = "list";
        String option1 = "x";
        String option2 = "y";

        List<StudyEvent> studyEvents = Lists.newArrayList(new StudyEvent(studyOID, Lists.newArrayList(formOID)));
        List<Form> forms = Lists.newArrayList(new Form(formOID, Lists.newArrayList(itemGroupOID)));
        List<ItemGroup> itemGroups = Lists.newArrayList(new ItemGroup(itemGroupOID, Lists.newArrayList(itemOID)));
        List<Item> items = Lists.newArrayList(new Item(itemOID, description, codeListOID));
        Map<Integer, String> codeListItems = Maps.newHashMap();
        codeListItems.put(1, option1);
        codeListItems.put(2, option2);

        List<CodeList> codeLists = Lists.newArrayList(new CodeList(codeListOID, codeListItems));
        XMLEcrfDatamodel datamodel = new XMLEcrfDatamodel(studyEvents, forms, itemGroups, items, codeLists);

        List<EcrfField> fields = XMLEcrfDatamodelToEcrfFields.convert(datamodel);

        assertEquals(1, fields.size());
        EcrfField field = fields.get(0);
        assertEquals("STUDY.FORM.ITEMGROUP.ITEM", field.name());
        assertEquals(description, field.description());

        Map<Integer, String> values = field.codeList();
        assertEquals(option1, values.get(1));
        assertEquals(option2, values.get(2));
    }
}