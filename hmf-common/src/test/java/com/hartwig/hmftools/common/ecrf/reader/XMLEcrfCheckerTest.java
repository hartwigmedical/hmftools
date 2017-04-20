package com.hartwig.hmftools.common.ecrf.reader;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.junit.Test;

public class XMLEcrfCheckerTest {
    private final String studyOID = "SE.Study";
    private final String formOID = "FRM.Form";
    private final String itemGroupOID = "GRP.ItemGroup";
    private final String itemOID = "GRP.Item";
    private final String description = "bla";
    private final String codeListOID = "list";
    private final String option1 = "x";
    private final String option2 = "y";

    @Test
    public void findsMissingCodeList() {
        final String codeListOIDRef = "CodeListOIDRef";

        final List<StudyEvent> studyEvents = Lists.newArrayList(new StudyEvent(studyOID, Lists.newArrayList(formOID)));
        final List<Form> forms = Lists.newArrayList(new Form(formOID, Lists.newArrayList(itemGroupOID)));
        final List<ItemGroup> itemGroups = Lists.newArrayList(
                new ItemGroup(itemGroupOID, Lists.newArrayList(itemOID)));
        final List<Item> items = Lists.newArrayList(new Item(itemOID, description, codeListOIDRef));
        final Map<Integer, String> codeListItems = Maps.newHashMap();
        codeListItems.put(1, option1);
        codeListItems.put(2, option2);

        final List<CodeList> codeLists = Lists.newArrayList(new CodeList(codeListOID, codeListItems));
        final XMLEcrfDatamodel datamodel = new XMLEcrfDatamodel(studyEvents, forms, itemGroups, items, codeLists);
        final List<String> fields = XMLEcrfChecker.checkReferences(datamodel);

        assertEquals(1, fields.size());
        final String missingItem = fields.get(0);
        assertEquals("Missing CodeList: " + codeListOIDRef, missingItem);
    }

    @Test
    public void findsMissingItem() {
        final String itemOIDRef = "ItemOIDRef";

        final List<StudyEvent> studyEvents = Lists.newArrayList(new StudyEvent(studyOID, Lists.newArrayList(formOID)));
        final List<Form> forms = Lists.newArrayList(new Form(formOID, Lists.newArrayList(itemGroupOID)));
        final List<ItemGroup> itemGroups = Lists.newArrayList(
                new ItemGroup(itemGroupOID, Lists.newArrayList(itemOIDRef)));
        final List<Item> items = Lists.newArrayList(new Item(itemOID, description, codeListOID));
        final Map<Integer, String> codeListItems = Maps.newHashMap();
        codeListItems.put(1, option1);
        codeListItems.put(2, option2);

        final List<CodeList> codeLists = Lists.newArrayList(new CodeList(codeListOID, codeListItems));
        final XMLEcrfDatamodel datamodel = new XMLEcrfDatamodel(studyEvents, forms, itemGroups, items, codeLists);
        final List<String> fields = XMLEcrfChecker.checkReferences(datamodel);

        assertEquals(1, fields.size());
        final String missingItem = fields.get(0);
        assertEquals("Missing Item: " + itemOIDRef, missingItem);
    }

    @Test
    public void findsMissingItemGroup() {
        final String itemGroupRef = "ItemGroupRef";

        final List<StudyEvent> studyEvents = Lists.newArrayList(new StudyEvent(studyOID, Lists.newArrayList(formOID)));
        final List<Form> forms = Lists.newArrayList(new Form(formOID, Lists.newArrayList(itemGroupRef)));
        final List<ItemGroup> itemGroups = Lists.newArrayList(
                new ItemGroup(itemGroupOID, Lists.newArrayList(itemOID)));
        final List<Item> items = Lists.newArrayList(new Item(itemOID, description, codeListOID));
        final Map<Integer, String> codeListItems = Maps.newHashMap();
        codeListItems.put(1, option1);
        codeListItems.put(2, option2);

        final List<CodeList> codeLists = Lists.newArrayList(new CodeList(codeListOID, codeListItems));
        final XMLEcrfDatamodel datamodel = new XMLEcrfDatamodel(studyEvents, forms, itemGroups, items, codeLists);
        final List<String> fields = XMLEcrfChecker.checkReferences(datamodel);

        assertEquals(1, fields.size());
        final String missingItem = fields.get(0);
        assertEquals("Missing ItemGroup: " + itemGroupRef, missingItem);
    }

    @Test
    public void findsMissingForm() {
        final String formRef = "formRef";

        final List<StudyEvent> studyEvents = Lists.newArrayList(new StudyEvent(studyOID, Lists.newArrayList(formRef)));
        final List<Form> forms = Lists.newArrayList(new Form(formOID, Lists.newArrayList(itemGroupOID)));
        final List<ItemGroup> itemGroups = Lists.newArrayList(
                new ItemGroup(itemGroupOID, Lists.newArrayList(itemOID)));
        final List<Item> items = Lists.newArrayList(new Item(itemOID, description, codeListOID));
        final Map<Integer, String> codeListItems = Maps.newHashMap();
        codeListItems.put(1, option1);
        codeListItems.put(2, option2);

        final List<CodeList> codeLists = Lists.newArrayList(new CodeList(codeListOID, codeListItems));
        final XMLEcrfDatamodel datamodel = new XMLEcrfDatamodel(studyEvents, forms, itemGroups, items, codeLists);
        final List<String> fields = XMLEcrfChecker.checkReferences(datamodel);

        assertEquals(1, fields.size());
        final String missingItem = fields.get(0);
        assertEquals("Missing Form: " + formRef, missingItem);
    }
}