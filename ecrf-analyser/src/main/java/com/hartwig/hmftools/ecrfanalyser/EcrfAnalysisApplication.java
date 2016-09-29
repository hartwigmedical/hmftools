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
        while (reader.hasNext() && count < 100000) {
            reader.next();
            if (reader.getEventType() == XMLEvent.START_ELEMENT) {
//                LOGGER.info(count + " Type=" + getEventTypeString(reader.getEventType()) + ", Name=" + reader.getName()
//                        + ", Count=" + reader.getAttributeCount());
                if (reader.getName().toString().equals("TranslatedText")) {
                    int x = 1;
                }
            }

            if (reader.hasText()) {
                LOGGER.info("Text: " + reader.getText());
            }
            count++;
        }

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


}
