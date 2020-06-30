package com.hartwig.hmftools.bachelor.types;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.BACH_LOGGER;

import java.nio.file.Path;

import javax.xml.XMLConstants;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.validation.Schema;
import javax.xml.validation.SchemaFactory;

import com.google.common.io.Resources;
import com.hartwig.hmftools.bachelor.datamodel.Program;

import org.jetbrains.annotations.Nullable;
import org.xml.sax.SAXException;

public class ConfigSchema
{
    private final Schema schema;

    private ConfigSchema() throws SAXException {
        schema = SchemaFactory.newInstance(XMLConstants.W3C_XML_SCHEMA_NS_URI).newSchema(Resources.getResource("bachelor.xsd"));
    }

    public static ConfigSchema make() throws SAXException {
        return new ConfigSchema();
    }

    @Nullable
    public Program processXML(final Path path)
    {
        BACH_LOGGER.info("Loading input file: {}", path);

        try
        {
            final JAXBContext context = JAXBContext.newInstance(Program.class);
            final Unmarshaller unmarshaller = context.createUnmarshaller();
            unmarshaller.setSchema(schema);
            return (Program) unmarshaller.unmarshal(path.toFile());
        }
        catch (final JAXBException e) {
            BACH_LOGGER.error("Failed to process: {}", path);
            BACH_LOGGER.error(e);
            return null;
        }
    }
}
