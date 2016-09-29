package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.List;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.events.XMLEvent;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.ecrfanalyser.EcrfAnalysisApplication;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfField;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class EcrfReader {

    private static final Logger LOGGER = LogManager.getLogger(EcrfReader.class);

    private static final String START_CLINICAL_DATA_TAG = "ClinicalData";
    private static final String START_ITEM_DEF_TAG = "ItemDef";
    private static final String ITEM_DEF_OID_ATTRIBUTE = "OID";
    private static final String ITEM_DEF_NAME_ATTRIBUTE = "Name";
    private static final String ITEM_DEF_CODE_LIST_REF = "CodeListRef";
    private static final String ITEM_DEF_CODE_LIST_OID = "CodeListOID";

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
            next(reader);
        }

        return new ODMContainer(itemDefs, codeLists);
    }

    @NotNull
    @VisibleForTesting
    static List<EcrfField> odmToEcrfFields(@NotNull ODMContainer container) {
        return Lists.newArrayList();

    }

    private static boolean isClinicalDataStart(@NotNull XMLStreamReader reader) {
        return reader.getEventType() == XMLEvent.START_ELEMENT && reader.getName().getLocalPart().equals(
                START_CLINICAL_DATA_TAG);
    }

    private static boolean isItemDefStart(@NotNull XMLStreamReader reader) {
        return reader.getEventType() == XMLEvent.START_ELEMENT && reader.getName().getLocalPart().equals(
                START_ITEM_DEF_TAG);
    }

    private static boolean isCodeListRefStart(@NotNull XMLStreamReader reader) {
        return reader.getEventType() == XMLEvent.START_ELEMENT && reader.getName().getLocalPart().equals(
                ITEM_DEF_CODE_LIST_REF);
    }

    @NotNull
    private static ItemDef extractItemDef(@NotNull XMLStreamReader reader) throws XMLStreamException {
        String OID = reader.getAttributeValue("", ITEM_DEF_OID_ATTRIBUTE);
        String name = reader.getAttributeValue("", ITEM_DEF_NAME_ATTRIBUTE);
        String codeListOID = null;
        while (reader.getEventType() != XMLEvent.END_ELEMENT) {
            next(reader);
            if (isCodeListRefStart(reader)) {
                codeListOID = reader.getAttributeValue("", ITEM_DEF_CODE_LIST_OID);
            }
        }

        return new ItemDef(OID, name, codeListOID);
    }

    private static void next(@NotNull XMLStreamReader reader) throws XMLStreamException {
        reader.next();
        String message = getEventTypeString(reader.getEventType());
        if (reader.getEventType() == XMLEvent.START_ELEMENT) {
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
