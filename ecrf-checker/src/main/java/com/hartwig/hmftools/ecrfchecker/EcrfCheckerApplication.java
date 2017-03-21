package com.hartwig.hmftools.ecrfchecker;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import com.hartwig.hmftools.common.ecrf.reader.XMLEcrfChecker;
import com.hartwig.hmftools.common.ecrf.reader.XMLEcrfDatamodel;
import com.hartwig.hmftools.common.ecrf.reader.XMLEcrfDatamodelReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class EcrfCheckerApplication {
    private static final Logger LOGGER = LogManager.getLogger(EcrfCheckerApplication.class);
    private static final String ECRF_XML_PATH = "ecrf";

    public static void main(final String... args) throws ParseException, IOException, XMLStreamException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        final String ecrfXmlPath = cmd.getOptionValue(ECRF_XML_PATH);

        if (ecrfXmlPath == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Ecrf-Checker", options);
            System.exit(1);
        }
        new EcrfCheckerApplication(ecrfXmlPath).run();
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(ECRF_XML_PATH, true, "The path to the ecrf xml file.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private final String ecrfXmlPath;

    private EcrfCheckerApplication(@NotNull final String ecrfXmlPath) {
        this.ecrfXmlPath = ecrfXmlPath;
    }

    private void run() throws IOException, XMLStreamException {
        final XMLInputFactory factory = XMLInputFactory.newInstance();
        final XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(ecrfXmlPath));
        final XMLEcrfDatamodel datamodel = XMLEcrfDatamodelReader.readXMLDatamodel(reader);
        final List<String> missingItems = XMLEcrfChecker.checkReferences(datamodel);
        for (String missingItem : missingItems) {
            LOGGER.info(missingItem);
        }
    }
}
