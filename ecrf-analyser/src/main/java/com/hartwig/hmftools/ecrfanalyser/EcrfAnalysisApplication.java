package com.hartwig.hmftools.ecrfanalyser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class EcrfAnalysisApplication {

    private static final Logger LOGGER = LogManager.getLogger(EcrfAnalysisApplication.class);

    public static void main(final String... args) {
        LOGGER.info("Hello World");
    }

    void run(@NotNull String path) throws FileNotFoundException, XMLStreamException {
        XMLInputFactory factory = XMLInputFactory.newInstance();
        XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(path));
        reader.close();
    }
}
