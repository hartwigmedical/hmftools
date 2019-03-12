package com.hartwig.hmftools.viccknowledgebase;

import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class jsonConverterApplication {
    private static final Logger LOGGER = LogManager.getLogger(jsonConverterApplication.class);

    private static final String ALL_JSON_FILE = "all.json";
    private static final String BRCA_JSON_FILE = "brca.json";
    private static final String CGI_JSON_FILE = "cgi.json";
    private static final String CIVIC_JSON_FILE = "civic.json";
    private static final String JAX_JSON_FILE = "jax.json";
    private static final String JAXTRIALS_JSON_FILE = "jax_trials.json";
    private static final String MOLECULARMATCH_JSON_FILE = "molecularmatch.json";
    private static final String MOLECULARMATCHTRIALS_JSON_FILE = "molecularmatch_trials.json";
    private static final String ONCOKB_JSON_FILE = "oncokb.json";
    private static final String PMKB_JSON_FILE = "pmkb.json";
    private static final String SAGE_JSON_FILE = "sage.json";

    public static void main(final String... args) throws ParseException, IOException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(ALL_JSON_FILE, true, "Path to json file of all knowledgebase together.");
        options.addOption(BRCA_JSON_FILE, true, "Path to the json file of brca knowledgebase.");
        options.addOption(CGI_JSON_FILE, true, "Path to the json file of cgi knowledgebase.");
        options.addOption(CIVIC_JSON_FILE, true, "Path to the json file of civic knowledgebase.");
        options.addOption(JAX_JSON_FILE, true, "Path to the json file of jax knowledgebase.");
        options.addOption(JAXTRIALS_JSON_FILE, true, "Path to the json file of jax trials knowledgebase.");
        options.addOption(MOLECULARMATCH_JSON_FILE, true, "Path to the json file of molecular match knowledgebase.");
        options.addOption(MOLECULARMATCHTRIALS_JSON_FILE, true, "Path to the json file of molecular match trials knowledgebase.");
        options.addOption(ONCOKB_JSON_FILE, true, "Path to the json file of oncokb knowledgebase.");
        options.addOption(PMKB_JSON_FILE, true, "Path to the json file of pmkb knowledgebase.");
        options.addOption(SAGE_JSON_FILE, true, "Path to the json file of sage knowledgebase.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static void printUsageAndExit(@NotNull final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("vicc-knowledgebase-importer", options);
        System.exit(1);
    }

}
