package com.hartwig.hmftools.serve;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.io.IclusionTrialFile;
import com.hartwig.hmftools.serve.sources.compassionateuse.CompassionateUseProgram;
import com.hartwig.hmftools.serve.sources.compassionateuse.CompassionateUseProgramFile;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.reader.ViccJsonReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ServeGeneratorApplication {

    private static final Logger LOGGER = LogManager.getLogger(ServeGeneratorApplication.class);
    private static final Set<String> PERMITTED_REF_GENOME_VERSIONS = Sets.newHashSet("hg19");

    private static final String VICC_JSON = "vicc_json";
    private static final String ICLUSION_TRIAL_TSV = "iclusion_trial_tsv";
    private static final String COMPASSIONATE_USE_PROGRAM_TSV = "compassionate_use_program_tsv";

    private static final String REF_GENOME_VERSION = "ref_genome_version";
    private static final String REF_GENOME_FASTA_FILE = "ref_genome_fasta_file";

    private static final String OUTPUT_DIR = "output_dir";

    private static final String VERSION = ServeGeneratorApplication.class.getPackage().getImplementationVersion();

    public static void main(String[] args) throws ParseException, IOException {
        LOGGER.info("Running Knowledgebase Generator v{}", VERSION);

        Options options = createOptions();
        CommandLine cmd = createCommandLine(args, options);

        if (!validInputForKnowledgebaseGeneration(cmd)) {
            printUsageAndExit(options);
        }

        // These are just to test the reading of the files. Handling will happen later.
        readIclusionTrials(cmd.getOptionValue(ICLUSION_TRIAL_TSV));
        readCompassionateUsePrograms(cmd.getOptionValue(COMPASSIONATE_USE_PROGRAM_TSV));

        List<ViccEntry> viccEntries = readViccEntries(cmd.getOptionValue(VICC_JSON));

        //TODO: generate HMF KB

        // Currently only support hg19.
        String refVersionString = cmd.getOptionValue(REF_GENOME_VERSION);
        assert PERMITTED_REF_GENOME_VERSIONS.contains(refVersionString);
        assert refVersionString.equals("hg19");
        RefGenomeVersion refGenomeVersion = RefGenomeVersion.HG19;

    }

    private static void readIclusionTrials(@NotNull String iClusionTrialTsv) throws IOException {
        LOGGER.info("Reading iClusion trials from {}", iClusionTrialTsv);
        List<IclusionTrial> trials = IclusionTrialFile.read(iClusionTrialTsv);
        LOGGER.info(" Read {} iClusion trials", trials.size());
    }

    private static void readCompassionateUsePrograms(@NotNull String compassionateUseProgramTsv) throws IOException {
        LOGGER.info("Reading compassionate use programs from {}", compassionateUseProgramTsv);
        List<CompassionateUseProgram> compassionateUsePrograms = CompassionateUseProgramFile.read(compassionateUseProgramTsv);
        LOGGER.info(" Read {} compassionate use programs", compassionateUsePrograms.size());
    }

    @NotNull
    private static List<ViccEntry> readViccEntries(@NotNull String viccJson) throws IOException {
        LOGGER.info("Reading VICC entries from {}", viccJson);
        List<ViccEntry> viccEntries = ViccJsonReader.buildProductionReader().readAll(viccJson);
        LOGGER.info(" Read {} VICC entries", viccEntries.size());
        return viccEntries;
    }

    private static boolean validInputForKnowledgebaseGeneration(@NotNull CommandLine cmd) {
        return fileExists(cmd, ICLUSION_TRIAL_TSV) && fileExists(cmd, VICC_JSON) && fileExists(cmd, COMPASSIONATE_USE_PROGRAM_TSV)
                && paramExists(cmd, REF_GENOME_VERSION) && valueIsPermitted(cmd, REF_GENOME_VERSION, PERMITTED_REF_GENOME_VERSIONS)
                && fileExists(cmd, REF_GENOME_FASTA_FILE) && dirExists(cmd, OUTPUT_DIR);
    }

    private static boolean fileExists(@NotNull CommandLine cmd, @NotNull String param) {
        if (paramExists(cmd, param) && !Files.exists(new File(cmd.getOptionValue(param)).toPath())) {
            LOGGER.warn("{} does not exist while '{}' has to be an existing path", cmd.getOptionValue(param), param);
            return false;
        }

        return true;
    }

    private static boolean dirExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);

        if (value == null || !pathExists(value) || !pathIsDirectory(value)) {
            LOGGER.warn(param + " has to be an existing directory: " + value);
            return false;
        }

        return true;
    }

    private static boolean pathExists(@NotNull String path) {
        return Files.exists(new File(path).toPath());
    }

    private static boolean pathIsDirectory(@NotNull String path) {
        return Files.isDirectory(new File(path).toPath());
    }

    private static boolean paramExists(@NotNull CommandLine cmd, @NotNull String param) {
        if (!cmd.hasOption(param)) {
            LOGGER.warn("Param '{}' has to be provided", param);
            return false;
        }

        return true;
    }

    private static boolean valueIsPermitted(@NotNull CommandLine cmd, @NotNull String param, @NotNull Set<String> permittedValues) {
        assert paramExists(cmd, param);
        String value = cmd.getOptionValue(param);

        if (!permittedValues.contains(value)) {
            LOGGER.warn("Value '{}' is not permitted for '{}'", value, param);
            return false;
        }

        return true;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        return new DefaultParser().parse(options, args);
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(VICC_JSON, true, "VICC JSON knowledgebase");
        options.addOption(ICLUSION_TRIAL_TSV, true, "iClusion input trial tsv");
        options.addOption(COMPASSIONATE_USE_PROGRAM_TSV, true, "Compassionate use program input tsv");

        options.addOption(REF_GENOME_VERSION, true, "Ref version. Should be 'hgxx'");
        options.addOption(REF_GENOME_FASTA_FILE, true, "Path to the ref genome fasta file");

        options.addOption(OUTPUT_DIR, true, "Path to the output dir of the files");

        return options;
    }

    private static void printUsageAndExit(@NotNull Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Knowledgebase-Generator", options);
        System.exit(1);
    }
}
