package com.hartwig.hmftools.bachelor;

import java.nio.file.Path;

import javax.xml.XMLConstants;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.validation.Schema;
import javax.xml.validation.SchemaFactory;

import com.google.common.io.Resources;

import nl.hartwigmedicalfoundation.bachelor.Program;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;
import org.xml.sax.SAXException;

public class BachelorSchema {

    private static final Logger LOGGER = LogManager.getLogger(BachelorSchema.class);
    private final Schema schema;

    private BachelorSchema() throws SAXException {
        schema = SchemaFactory.newInstance(XMLConstants.W3C_XML_SCHEMA_NS_URI).newSchema(Resources.getResource("bachelor.xsd"));
    }

    public static BachelorSchema make() throws SAXException {
        return new BachelorSchema();
    }

    @Nullable
    public Program processXML(final Path path) {
        LOGGER.info("loading file: {}", path);
        try {
            final JAXBContext context = JAXBContext.newInstance(Program.class);
            final Unmarshaller unmarshaller = context.createUnmarshaller();
            unmarshaller.setSchema(schema);
            return (Program) unmarshaller.unmarshal(path.toFile());
        } catch (final JAXBException e) {
            LOGGER.error("failed to process: {}", path);
            LOGGER.error(e);
            return null;
        }
    }
}
