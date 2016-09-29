package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.events.XMLEvent;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfField;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class EcrfReader {

    private static final Logger LOGGER = LogManager.getLogger(EcrfReader.class);

    private static final String CLINICAL_DATA_TAG = "ClinicalData";

    private static final String ITEM_DEF_TAG = "ItemDef";
    private static final String ITEM_DEF_OID_ATTRIBUTE = "OID";
    private static final String ITEM_DEF_NAME_ATTRIBUTE = "Name";
    private static final String ITEM_DEF_CODE_LIST_REF = "CodeListRef";
    private static final String ITEM_DEF_CODE_LIST_OID = "CodeListOID";

    private static final String CODE_LIST_TAG = "CodeList";
    private static final String CODE_LIST_OID_ATTRIBUTE = "OID";
    private static final String CODE_LIST_ITEM = "TranslatedText";

    private EcrfReader() {
    }

    @NotNull
    static List<EcrfField> extractFields(@NotNull XMLStreamReader reader) throws XMLStreamException {
        ODMContainer container = extractODM(reader);
        return odmToEcrfFields(container);
    }

    @NotNull
    @VisibleForTesting
    static ODMContainer extractODM(@NotNull XMLStreamReader reader) throws XMLStreamException {
        List<ItemDef> itemDefs = Lists.newArrayList();
        List<CodeList> codeLists = Lists.newArrayList();
        while (reader.hasNext() && !isClinicalDataStart(reader)) {
            if (isItemDefStart(reader)) {
                itemDefs.add(extractItemDef(reader));
            }
            if (isCodeListStart(reader)) {
                codeLists.add(extractCodeList(reader));
            }
            next(reader);
        }

        return new ODMContainer(itemDefs, codeLists);
    }

    @NotNull
    @VisibleForTesting
    static List<EcrfField> odmToEcrfFields(@NotNull ODMContainer container) {
        List<EcrfField> fields = Lists.newArrayListWithCapacity(container.itemDefs().size());
        for (ItemDef item : container.itemDefs()) {
            Map<Integer, String> values = Maps.newHashMap();
            if (item.codeListOID() != null) {
                values = findValuesForCodeList(container.codeLists(), item.codeListOID());
            }
            String category = ItemDefToEcrfField.category(item);
            String fieldName = ItemDefToEcrfField.fieldName(item);
            String description = ItemDefToEcrfField.description(item);

            fields.add(new EcrfField(category, fieldName, description, values));
        }
        return fields;
    }

    @NotNull
    private static Map<Integer,String> findValuesForCodeList(final List<CodeList> codeLists, final String codeListOID) {
        for (CodeList codeList : codeLists) {
            if (codeList.OID().equals(codeListOID)) {
                return codeList.values();
            }
        }

        throw new IllegalStateException("Could not find code list : " + codeListOID);
    }

    @NotNull
    private static ItemDef extractItemDef(@NotNull XMLStreamReader reader) throws XMLStreamException {
        String OID = reader.getAttributeValue("", ITEM_DEF_OID_ATTRIBUTE);
        String name = reader.getAttributeValue("", ITEM_DEF_NAME_ATTRIBUTE);
        String codeListOID = null;
        while (!isItemDefEnd(reader)) {
            next(reader);
            if (isCodeListRefStart(reader)) {
                codeListOID = reader.getAttributeValue("", ITEM_DEF_CODE_LIST_OID);
            }
        }

        return new ItemDef(OID, name, codeListOID);
    }

    @NotNull
    private static CodeList extractCodeList(@NotNull XMLStreamReader reader) throws XMLStreamException {
        String OID = reader.getAttributeValue("", CODE_LIST_OID_ATTRIBUTE);
        List<String> codeListItems = Lists.newArrayList();
        while (!isCodeListEnd(reader)) {
            next(reader);
            if (isCodeListItem(reader)) {
                // KODU: The item content for the code lists is held in text in the next event!
                next(reader);
                codeListItems.add(reader.getText());
            }
        }
        return new CodeList(OID, CodeListFactory.extractValuesFromStrings(codeListItems));
    }

    private static boolean isClinicalDataStart(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, CLINICAL_DATA_TAG);
    }

    private static boolean isItemDefStart(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, ITEM_DEF_TAG);
    }

    private static boolean isItemDefEnd(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.END_ELEMENT, ITEM_DEF_TAG);
    }

    private static boolean isCodeListRefStart(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, ITEM_DEF_CODE_LIST_REF);
    }

    private static boolean isCodeListStart(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, CODE_LIST_TAG);
    }

    private static boolean isCodeListEnd(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.END_ELEMENT, CODE_LIST_TAG);
    }

    private static boolean isCodeListItem(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, CODE_LIST_ITEM);
    }


    private static boolean isOfTypeWithName(@NotNull XMLStreamReader reader, int event,
            @NotNull String name) {
        return reader.getEventType() == event && reader.getName().getLocalPart().equals(name);
    }

    private static void next(@NotNull XMLStreamReader reader) throws XMLStreamException {
        reader.next();
        String message = getEventTypeString(reader.getEventType());
        if (reader.getEventType() == XMLEvent.START_ELEMENT || reader.getEventType() == XMLEvent.END_ELEMENT) {
            message += ": " + reader.getName().toString();
        }
        LOGGER.info(message);
    }

    @NotNull
    private static String getEventTypeString(int eventType) {
        switch (eventType) {
            case XMLEvent.START_ELEMENT:
                return "START_ELEMENT";
            case XMLEvent.END_ELEMENT:
                return "END_ELEMENT";
            case XMLEvent.PROCESSING_INSTRUCTION:
                return "PROCESSING_INSTRUCTION";
            case XMLEvent.CHARACTERS:
                return "CHARACTERS";
            case XMLEvent.COMMENT:
                return "COMMENT";
            case XMLEvent.START_DOCUMENT:
                return "START_DOCUMENT";
            case XMLEvent.END_DOCUMENT:
                return "END_DOCUMENT";
            case XMLEvent.ENTITY_REFERENCE:
                return "ENTITY_REFERENCE";
            case XMLEvent.ATTRIBUTE:
                return "ATTRIBUTE";
            case XMLEvent.DTD:
                return "DTD";
            case XMLEvent.CDATA:
                return "CDATA";
            case XMLEvent.SPACE:
                return "SPACE";
        }
        return "UNKNOWN_EVENT_TYPE , " + eventType;
    }
}
