package com.hartwig.hmftools.common.ecrf.reader;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfField;

import org.junit.Test;

public class XMLEcrfDatamodelToEcrfFieldsTest {

    @Test
    public void canConvertXMLObjectContainerToEcrfFields() {
        final String studyOID = "SE.Study";
        final String formOID = "FRM.Form";
        final String itemGroupOID = "GRP.ItemGroup";
        final String itemOID = "GRP.Item";
        final String description = "bla";

        final String codeListOID = "list";
        final String option1 = "x";
        final String option2 = "y";

        final List<StudyEvent> studyEvents =
                Lists.newArrayList(new ImmutableStudyEvent(studyOID, description, Lists.newArrayList(formOID)));
        final List<Form> forms = Lists.newArrayList(new ImmutableForm(formOID, description, Lists.newArrayList(itemGroupOID)));
        final List<ItemGroup> itemGroups =
                Lists.newArrayList(new ImmutableItemGroup(itemGroupOID, description, Lists.newArrayList(itemOID)));
        final List<Item> items = Lists.newArrayList(new ImmutableItem(itemOID, description, codeListOID));
        final Map<Integer, String> codeListItems = Maps.newHashMap();
        codeListItems.put(1, option1);
        codeListItems.put(2, option2);

        final List<CodeList> codeLists = Lists.newArrayList(new ImmutableCodeList(codeListOID, description, codeListItems));
        final XMLEcrfDatamodel datamodel = new ImmutableXMLEcrfDatamodel(studyEvents, forms, itemGroups, items, codeLists);

        final List<EcrfField> fields = XMLEcrfDatamodelToEcrfFields.convert(datamodel);

        assertEquals(1, fields.size());
        final EcrfField field = fields.get(0);
        assertEquals("STUDY.FORM.ITEMGROUP.ITEM", field.name());
        assertEquals(description, field.description());

        final Map<Integer, String> values = field.codeList();
        assertEquals(option1, values.get(1));
        assertEquals(option2, values.get(2));
    }
}