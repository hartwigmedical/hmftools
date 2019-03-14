package com.hartwig.hmftools.viccknowledgebase;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import com.hartwig.hmftools.common.vicc.ViccFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class JsonConverterApplication {
    private static final Logger LOGGER = LogManager.getLogger(JsonConverterApplication.class);

    private static final String ALL_JSON_FILE = "allFile";
    private static final String BRCA_JSON_FILE = "brcaFile";
    private static final String CGI_JSON_FILE = "cgiFile";
    private static final String CIVIC_JSON_FILE = "civicFile";
    private static final String JAX_JSON_FILE = "jaxFile";
    private static final String JAXTRIALS_JSON_FILE = "jaxTrialsFile";
    private static final String MOLECULARMATCH_JSON_FILE = "molecularMatchFile";
    private static final String MOLECULARMATCHTRIALS_JSON_FILE = "molecularMatchTrialsFile";
    private static final String ONCOKB_JSON_FILE = "oncokbFile";
    private static final String PMKB_JSON_FILE = "pmkbFile";
    private static final String SAGE_JSON_FILE = "sageFile";
    private static final String PATH_KNOWLEDGEBASE_FILES = "knowledgebase";

    public static void main(final String... args) throws ParseException, IOException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        if (!validInputForBaseReportData(cmd)) {
            printUsageAndExit(options);
        }

        LOGGER.info("Running json converter to csv");

        String knowledgebasePath = cmd.getOptionValue(PATH_KNOWLEDGEBASE_FILES);

        readingAllJsonFile(cmd.getOptionValue(ALL_JSON_FILE), knowledgebasePath);
        readingBrcaJsonFile(cmd.getOptionValue(BRCA_JSON_FILE), knowledgebasePath);
        readingCgiJsonFile(cmd.getOptionValue(CGI_JSON_FILE), knowledgebasePath);
        readingCivicJsonFile(cmd.getOptionValue(CIVIC_JSON_FILE), knowledgebasePath);
        readingJaxJsonFile(cmd.getOptionValue(JAX_JSON_FILE), knowledgebasePath);
        readingJaxTrialsJsonFile(cmd.getOptionValue(JAXTRIALS_JSON_FILE), knowledgebasePath);
        readingMolecularMatchJsonFile(cmd.getOptionValue(MOLECULARMATCH_JSON_FILE), knowledgebasePath);
        readingMolecularMatchTrialsJsonFile(cmd.getOptionValue(MOLECULARMATCHTRIALS_JSON_FILE), knowledgebasePath);
        readingOncokbJsonFile(cmd.getOptionValue(ONCOKB_JSON_FILE), knowledgebasePath);
        readingPmkbJsonFile(cmd.getOptionValue(PMKB_JSON_FILE), knowledgebasePath);
        readingSageJsonFile(cmd.getOptionValue(SAGE_JSON_FILE), knowledgebasePath);
    }

    private static void readingAllJsonFile(@NotNull String allJsonFile, @NotNull String knowledgebasePath) throws IOException {
        String allJsonFilePath = knowledgebasePath + File.separator + allJsonFile;
        ViccFactory.extractAllFile(allJsonFilePath);
    }

    private static void readingBrcaJsonFile(@NotNull String brcaJsonFile, @NotNull String knowledgebasePath) throws IOException {
        String brcaJsonFilePath = knowledgebasePath + File.separator + brcaJsonFile;
        ViccFactory.extractBRCAFile(brcaJsonFilePath);
    }

    private static void readingCgiJsonFile(@NotNull String cgiJsonFile, @NotNull String knowledgebasePath) throws IOException {
        String cgiJsonFilePath = knowledgebasePath + File.separator + cgiJsonFile;
        ViccFactory.extractCgiFile(cgiJsonFilePath);
    }

    private static void readingCivicJsonFile(@NotNull String civicJsonFile, @NotNull String knowledgebasePath) throws IOException {
        String civicJsonFilePath = knowledgebasePath + File.separator + civicJsonFile;
        ViccFactory.extractCivicFile(civicJsonFilePath);
    }

    private static void readingJaxJsonFile(@NotNull String jaxJsonFile, @NotNull String knowledgebasePath) throws IOException {
        String jaxJsonFilePath = knowledgebasePath + File.separator + jaxJsonFile;
        ViccFactory.extractJaxFile(jaxJsonFilePath);
    }

    private static void readingJaxTrialsJsonFile(@NotNull String jaxTrialsJsonFile, @NotNull String knowledgebasePath) throws IOException {
        String jaxTrialsJsonFilePath = knowledgebasePath + File.separator + jaxTrialsJsonFile;
        ViccFactory.extractJaxTrialsFile(jaxTrialsJsonFilePath);
    }

    private static void readingMolecularMatchJsonFile(@NotNull String molecularMatchJsonFile, @NotNull String knowledgebasePath)
            throws IOException {
        String molecularMatchJsonFilePath = knowledgebasePath + File.separator + molecularMatchJsonFile;
        ViccFactory.extractMolecularMatchFile(molecularMatchJsonFilePath);
    }

    private static void readingMolecularMatchTrialsJsonFile(@NotNull String molecularMatchTrialsJsonFile, @NotNull String knowledgebasePath)
            throws IOException {
        String molecularMatchTrialsJsonFilePath = knowledgebasePath + File.separator + molecularMatchTrialsJsonFile;
        ViccFactory.extractMolecularMatchTrailsFile(molecularMatchTrialsJsonFilePath);
    }

    private static void readingOncokbJsonFile(@NotNull String oncokbJsonFile, @NotNull String knowledgebasePath) throws IOException {
        String oncokbJsonFilePath = knowledgebasePath + File.separator + oncokbJsonFile;
        ViccFactory.extractOncokbFile(oncokbJsonFilePath);
    }

    private static void readingPmkbJsonFile(@NotNull String pmkbJsonFile, @NotNull String knowledgebasePath) throws IOException {
        String pmkbJsonFilePath = knowledgebasePath + File.separator + pmkbJsonFile;
        ViccFactory.extractPmkbFile(pmkbJsonFilePath);
    }

    private static void readingSageJsonFile(@NotNull String sageJsonFile, @NotNull String knowledgebasePath) throws IOException {
        String sageJsonFilePath = knowledgebasePath + File.separator + sageJsonFile;
        ViccFactory.extractSageFile(sageJsonFilePath);
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
        options.addOption(PATH_KNOWLEDGEBASE_FILES, true, "Complete path a directory holding the knowledgebase data");
        return options;
    }

    private static boolean validInputForBaseReportData(@NotNull final CommandLine cmd) {
        final String allJson = cmd.getOptionValue(ALL_JSON_FILE);
        final String brcaJson = cmd.getOptionValue(BRCA_JSON_FILE);
        final String cgiJson = cmd.getOptionValue(CGI_JSON_FILE);
        final String jaxJson = cmd.getOptionValue(JAX_JSON_FILE);
        final String jaxTrialsJson = cmd.getOptionValue(JAXTRIALS_JSON_FILE);
        final String molecularMatchJson = cmd.getOptionValue(MOLECULARMATCH_JSON_FILE);
        final String molecularMatchTrialsJson = cmd.getOptionValue(MOLECULARMATCHTRIALS_JSON_FILE);
        final String oncokbJson = cmd.getOptionValue(ONCOKB_JSON_FILE);
        final String pmkbJson = cmd.getOptionValue(PMKB_JSON_FILE);
        final String sageJson = cmd.getOptionValue(SAGE_JSON_FILE);
        final String pathKnowledgebaseFiles = cmd.getOptionValue(PATH_KNOWLEDGEBASE_FILES);

        if (allJson == null || !exists(allJson)) {
            LOGGER.warn(ALL_JSON_FILE + " has to be an existing file: " + allJson);
        } else if (brcaJson == null || !exists(brcaJson)) {
            LOGGER.warn(BRCA_JSON_FILE + " has to be an existing file: " + brcaJson);
        } else if (cgiJson == null || !exists(cgiJson)) {
            LOGGER.warn(CGI_JSON_FILE + " has to be an existing file: " + cgiJson);
        } else if (jaxJson == null || !exists(jaxJson)) {
            LOGGER.warn(JAX_JSON_FILE + " has to be an existing file: " + jaxJson);
        } else if (jaxTrialsJson == null || !exists(jaxTrialsJson)) {
            LOGGER.warn(JAXTRIALS_JSON_FILE + " has to be an existing file: " + jaxTrialsJson);
        } else if (molecularMatchJson == null || !exists(molecularMatchJson)) {
            LOGGER.warn(MOLECULARMATCH_JSON_FILE + " has to be an existing file: " + molecularMatchJson);
        } else if (molecularMatchTrialsJson == null || !exists(molecularMatchTrialsJson)) {
            LOGGER.warn(MOLECULARMATCHTRIALS_JSON_FILE + " has to be an existing file: " + molecularMatchTrialsJson);
        } else if (oncokbJson == null || !exists(oncokbJson)) {
            LOGGER.info(ONCOKB_JSON_FILE + " has to be an existing file: " + oncokbJson);
        } else if (pmkbJson == null || !exists(pmkbJson)) {
            LOGGER.warn(PMKB_JSON_FILE + " has to be an existing file: " + pmkbJson);
        } else if (sageJson == null || !exists(sageJson)) {
            LOGGER.warn(SAGE_JSON_FILE + " has to be an existing file: " + sageJson);
        } else if (pathKnowledgebaseFiles == null || !exists(pathKnowledgebaseFiles) || !isDirectory(pathKnowledgebaseFiles)) {
            LOGGER.warn(PATH_KNOWLEDGEBASE_FILES + " has to be an existing directory: " + pathKnowledgebaseFiles);
        } else {
            return true;
        }
        return false;
    }

    private static boolean exists(@NotNull final String path) {
        return Files.exists(new File(path).toPath());
    }

    private static boolean isDirectory(@NotNull final String path) {
        return Files.isDirectory(new File(path).toPath());
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
