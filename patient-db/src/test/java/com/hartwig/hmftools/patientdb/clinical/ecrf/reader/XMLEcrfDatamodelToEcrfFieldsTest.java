package com.hartwig.hmftools.patientdb.clinical.ecrf.reader;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfDatamodelField;

import org.junit.Test;

public class XMLEcrfDatamodelToEcrfFieldsTest {

    @Test
    public void canConvertXMLObjectContainerToEcrfFields() {
        String studyOID = "SE.Study";
        String formOID = "FRM.Form";
        String itemGroupOID = "GRP.ItemGroup";
        String itemOID = "GRP.Item";
        String description = "bla";

        String codeListOID = "list";
        String option1 = "x";
        String option2 = "y";

        List<StudyEvent> studyEvents = Lists.newArrayList(new ImmutableStudyEvent(studyOID, description, Lists.newArrayList(formOID)));
        List<Form> forms = Lists.newArrayList(new ImmutableForm(formOID, description, Lists.newArrayList(itemGroupOID)));
        List<ItemGroup> itemGroups = Lists.newArrayList(new ImmutableItemGroup(itemGroupOID, description, Lists.newArrayList(itemOID)));
        List<Item> items = Lists.newArrayList(new ImmutableItem(itemOID, description, codeListOID));
        Map<Integer, String> codeListItems = Maps.newHashMap();
        codeListItems.put(1, option1);
        codeListItems.put(2, option2);

        List<CodeList> codeLists = Lists.newArrayList(new ImmutableCodeList(codeListOID, description, codeListItems));
        XMLEcrfDatamodel datamodel = XMLEcrfDatamodel.of(studyEvents, forms, itemGroups, items, codeLists);

        List<EcrfDatamodelField> fields = XMLEcrfDatamodelToEcrfFields.convert(datamodel);

        assertEquals(1, fields.size());
        EcrfDatamodelField field = fields.get(0);
        assertEquals("STUDY.FORM.ITEMGROUP.ITEM", field.name());
        assertEquals(description, field.description());

        Map<Integer, String> values = field.codeList();
        assertEquals(option1, values.get(1));
        assertEquals(option2, values.get(2));
    }
}
