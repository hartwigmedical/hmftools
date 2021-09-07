package com.hartwig.hmftools.patientdb.clinical.ecrf.reader;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.events.XMLEvent;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

abstract class EcrfReader {

    private static final Logger LOGGER = LogManager.getLogger(EcrfReader.class);
    private static final boolean LOG_EVENTS = false;

    private static final String CLINICAL_DATA_TAG = "ClinicalData";

    static boolean isClinicalDataStart(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, CLINICAL_DATA_TAG);
    }

    static boolean isClinicalDataEnd(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.END_ELEMENT, CLINICAL_DATA_TAG);
    }

    static boolean isOfTypeWithName(@NotNull XMLStreamReader reader, int event, @NotNull String name) {
        return reader.getEventType() == event && reader.getName().getLocalPart().equals(name);
    }

    static void next(@NotNull XMLStreamReader reader) throws XMLStreamException {
        reader.next();
        if (LOG_EVENTS) {
            String message = getEventTypeString(reader.getEventType());
            if (reader.getEventType() == XMLEvent.START_ELEMENT || reader.getEventType() == XMLEvent.END_ELEMENT) {
                message += ": " + reader.getName().toString();
            }
            LOGGER.info(message);
        }
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
