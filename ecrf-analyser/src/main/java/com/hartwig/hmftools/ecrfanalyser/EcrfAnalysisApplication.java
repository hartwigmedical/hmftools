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

import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfField;
import com.hartwig.hmftools.ecrfanalyser.reader.EcrfReader;

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

        List<EcrfField> datamodel = EcrfReader.extractFields(reader);
        int x = 1;

    }

}
