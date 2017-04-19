package com.hartwig.hmftools.batchvariantanalyser;

import java.io.IOException;

import javax.xml.stream.XMLStreamException;

import com.hartwig.hmftools.common.exception.HartwigException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class BatchVariantAnalyserApplication {

    private static final Logger LOGGER = LogManager.getLogger(BatchVariantAnalyserApplication.class);

    private static final String INPUT_DIRECTORY = "input";
    private static final String OUTPUT_DIRECTORY = "output";

    public static void main(final String... args)
            throws ParseException, IOException, XMLStreamException, HartwigException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        LOGGER.info(cmd.getArgs());
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(INPUT_DIRECTORY, true, "The full path towards the input directory");
        options.addOption(OUTPUT_DIRECTORY, true, "The full path towards the output directory");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
