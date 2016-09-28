package com.hartwig.hmftools.ecrfanalyser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.List;

import javax.xml.stream.XMLEventFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.events.XMLEvent;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class EcrfAnalysisApplication {

    private static final Logger LOGGER = LogManager.getLogger(EcrfAnalysisApplication.class);

    public static void main(final String... args) {
        LOGGER.info("Hello World");
    }

    @NotNull
    private final String ecrfXmlPath;
    @NotNull
    private final String csvOutPath;

    EcrfAnalysisApplication(@NotNull final String ecrfXmlPath, @NotNull final String csvOutPath) {
        this.ecrfXmlPath = ecrfXmlPath;
        this.csvOutPath = csvOutPath;
    }

    void generateCsv(@NotNull List<String> patients, @NotNull List<String> fields)
            throws FileNotFoundException, XMLStreamException {
        XMLInputFactory factory = XMLInputFactory.newInstance();
        XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(ecrfXmlPath));
        int count = 0;
        while (reader.hasNext() && count < 10) {
            reader.next();
            LOGGER.info("ElementType " + count + ": " + getEventTypeString(reader.getEventType()));
            if (reader.getEventType() != XMLEvent.CHARACTERS) {
                LOGGER.info("ElementName " + count + ": " + reader.getName());
            }

            if (reader.getEventType() == XMLEvent.START_ELEMENT) {
                if (reader.hasText()) {
                    LOGGER.info(reader.getElementText());
                }
                LOGGER.info(reader.getAttributeCount());
            }

            reader.next();
            count++;
        }
        LOGGER.info(count);

//        int namespaces = reader.getNamespaceCount();
//        for (int i = 0; i < namespaces; i++) {
//            LOGGER.info("NamespacePrefix " + i + ": " + reader.getNamespacePrefix(i));
//            LOGGER.info("NamespaceURI " + i + ": " + reader.getNamespaceURI(i));
//        }
        reader.close();

        //        XMLEventReader eventReader = factory.createXMLEventReader(new FileInputStream(ecrfXmlPath));
        //        while (eventReader.hasNext()) {
        //            XMLEvent e = eventReader.nextEvent();
        //            LOGGER.info(e.toString());
        //        }
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
