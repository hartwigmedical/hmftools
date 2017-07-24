package com.hartwig.hmftools.common.ecrf.reader;

import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.events.XMLEvent;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public final class XMLEcrfDatamodelReader extends EcrfReader {

    private static final String STUDY_EVENT_TAG = "StudyEventDef";
    private static final String STUDY_EVENT_OID_ATTRIBUTE = "OID";
    private static final String STUDY_EVENT_FORM_REF = "FormRef";
    private static final String STUDY_EVENT_FORM_OID = "FormOID";

    private static final String FORM_TAG = "FormDef";
    private static final String FORM_OID_ATTRIBUTE = "OID";
    private static final String FORM_ITEM_GROUP_REF = "ItemGroupRef";
    private static final String FORM_ITEM_GROUP_OID = "ItemGroupOID";
    // KODU: This item group is referenced from forms but lacks definition!
    private static final String FORM_ITEM_GROUP_OID_IGNORE = "GRP.AuditData";

    private static final String ITEM_GROUP_TAG = "ItemGroupDef";
    private static final String ITEM_GROUP_OID_ATTRIBUTE = "OID";
    private static final String ITEM_GROUP_ITEM_REF = "ItemRef";
    private static final String ITEM_GROUP_ITEM_OID = "ItemOID";

    private static final String ITEM_TAG = "ItemDef";
    private static final String ITEM_OID_ATTRIBUTE = "OID";
    private static final String ITEM_NAME_ATTRIBUTE = "Name";
    private static final String ITEM_CODE_LIST_REF = "CodeListRef";
    private static final String ITEM_CODE_LIST_OID = "CodeListOID";

    private static final String CODE_LIST_TAG = "CodeList";
    private static final String CODE_LIST_OID_ATTRIBUTE = "OID";
    private static final String CODE_LIST_ITEM_TAG = "CodeListItem";
    private static final String CODE_LIST_ITEM_VALUE_ATTRIBUTE = "CodedValue";
    private static final String CODE_LIST_ITEM_CONTENT = "TranslatedText";

    private XMLEcrfDatamodelReader() {
    }

    @NotNull
    public static XMLEcrfDatamodel readXMLDatamodel(@NotNull final XMLStreamReader reader) throws XMLStreamException {
        final List<StudyEvent> studyEvents = Lists.newArrayList();
        final List<Form> forms = Lists.newArrayList();
        final List<ItemGroup> itemGroups = Lists.newArrayList();
        final List<Item> items = Lists.newArrayList();
        final List<CodeList> codeLists = Lists.newArrayList();
        while (reader.hasNext() && !isClinicalDataStart(reader)) {
            if (isStudyEventStart(reader)) {
                studyEvents.add(extractStudyEvent(reader));
            } else if (isFormStart(reader)) {
                forms.add(extractForm(reader));
            } else if (isItemGroupStart(reader)) {
                itemGroups.add(extractItemGroup(reader));
            } else if (isItemStart(reader)) {
                items.add(extractItem(reader));
            } else if (isCodeListStart(reader)) {
                codeLists.add(extractCodeList(reader));
            }
            next(reader);
        }

        return XMLEcrfDatamodel.of(studyEvents, forms, itemGroups, items, codeLists);
    }

    @NotNull
    private static StudyEvent extractStudyEvent(@NotNull final XMLStreamReader reader) throws XMLStreamException {
        final String OID = reader.getAttributeValue("", STUDY_EVENT_OID_ATTRIBUTE);
        final String name = reader.getAttributeValue("", ITEM_NAME_ATTRIBUTE);
        final List<String> subOIDs = Lists.newArrayList();
        while (!isStudyEventEnd(reader)) {
            next(reader);
            if (isFormRefStart(reader)) {
                subOIDs.add(reader.getAttributeValue("", STUDY_EVENT_FORM_OID));
            }
        }

        return new ImmutableStudyEvent(OID, name, subOIDs);
    }

    @NotNull
    private static Form extractForm(@NotNull final XMLStreamReader reader) throws XMLStreamException {
        final String OID = reader.getAttributeValue("", FORM_OID_ATTRIBUTE);
        final String name = reader.getAttributeValue("", ITEM_NAME_ATTRIBUTE);
        final List<String> subOIDs = Lists.newArrayList();
        while (!isFormEnd(reader)) {
            next(reader);
            if (isItemGroupRefStart(reader)) {
                String subOID = reader.getAttributeValue("", FORM_ITEM_GROUP_OID);
                if (!subOID.equals(FORM_ITEM_GROUP_OID_IGNORE)) {
                    subOIDs.add(subOID);
                }
            }
        }

        return new ImmutableForm(OID, name, subOIDs);
    }

    @NotNull
    private static ItemGroup extractItemGroup(@NotNull final XMLStreamReader reader) throws XMLStreamException {
        final String OID = reader.getAttributeValue("", ITEM_GROUP_OID_ATTRIBUTE);
        final String name = reader.getAttributeValue("", ITEM_NAME_ATTRIBUTE);
        final List<String> subOIDs = Lists.newArrayList();
        while (!isItemGroupEnd(reader)) {
            next(reader);
            if (isItemRefStart(reader)) {
                subOIDs.add(reader.getAttributeValue("", ITEM_GROUP_ITEM_OID));
            }
        }

        return new ImmutableItemGroup(OID, name, subOIDs);
    }

    @NotNull
    private static Item extractItem(@NotNull final XMLStreamReader reader) throws XMLStreamException {
        final String OID = reader.getAttributeValue("", ITEM_OID_ATTRIBUTE);
        final String name = reader.getAttributeValue("", ITEM_NAME_ATTRIBUTE);
        String codeListOID = null;
        while (!isItemEnd(reader)) {
            next(reader);
            if (isCodeListRefStart(reader)) {
                codeListOID = reader.getAttributeValue("", ITEM_CODE_LIST_OID);
            }
        }

        return new ImmutableItem(OID, name, codeListOID);
    }

    @NotNull
    private static CodeList extractCodeList(@NotNull final XMLStreamReader reader) throws XMLStreamException {
        final String OID = reader.getAttributeValue("", CODE_LIST_OID_ATTRIBUTE);
        final String name = reader.getAttributeValue("", ITEM_NAME_ATTRIBUTE);
        final Map<Integer, String> codeListItems = Maps.newHashMap();
        while (!isCodeListEnd(reader)) {
            next(reader);
            if (isCodeListItem(reader)) {
                int index = Integer.valueOf(reader.getAttributeValue("", CODE_LIST_ITEM_VALUE_ATTRIBUTE));
                while (!isCodeListItemContent(reader)) {
                    next(reader);
                }
                // KODU: The item content for the code lists is held in text in the next event!
                next(reader);
                codeListItems.put(index, CodeListFactory.fromText(reader.getText()));
            }
        }
        return new ImmutableCodeList(OID, name, codeListItems);
    }

    private static boolean isStudyEventStart(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, STUDY_EVENT_TAG);
    }

    private static boolean isFormRefStart(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, STUDY_EVENT_FORM_REF);
    }

    private static boolean isStudyEventEnd(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.END_ELEMENT, STUDY_EVENT_TAG);
    }

    private static boolean isFormStart(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, FORM_TAG);
    }

    private static boolean isItemGroupRefStart(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, FORM_ITEM_GROUP_REF);
    }

    private static boolean isFormEnd(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.END_ELEMENT, FORM_TAG);
    }

    private static boolean isItemGroupStart(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, ITEM_GROUP_TAG);
    }

    private static boolean isItemRefStart(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, ITEM_GROUP_ITEM_REF);
    }

    private static boolean isItemGroupEnd(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.END_ELEMENT, ITEM_GROUP_TAG);
    }

    private static boolean isItemStart(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, ITEM_TAG);
    }

    private static boolean isCodeListRefStart(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, ITEM_CODE_LIST_REF);
    }

    private static boolean isItemEnd(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.END_ELEMENT, ITEM_TAG);
    }

    private static boolean isCodeListStart(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, CODE_LIST_TAG);
    }

    private static boolean isCodeListEnd(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.END_ELEMENT, CODE_LIST_TAG);
    }

    private static boolean isCodeListItem(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, CODE_LIST_ITEM_TAG);
    }

    private static boolean isCodeListItemContent(@NotNull final XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, CODE_LIST_ITEM_CONTENT);
    }
}
