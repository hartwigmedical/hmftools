package com.hartwig.hmftools.knowledgebasegenerator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.knowledgebasegenerator.CNV.GeneratingCNV;
import com.hartwig.hmftools.knowledgebasegenerator.transvar.RefVersion;
import com.hartwig.hmftools.knowledgebasegenerator.transvar.Transvar;
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

public class KnowledgebaseGeneratorApplication {

    private static final Logger LOGGER = LogManager.getLogger(KnowledgebaseGeneratorApplication.class);

    private static final String VICC_JSON = "vicc_json";
    private static final String ICLUSION_TRIAL_TSV = "iclusion_trial_tsv";
    private static final String COMPASSIONATE_USE_PROGRAMS_TSV = "compassionate_use_programs_tsv";

    private static final String VERSION = KnowledgebaseGeneratorApplication.class.getPackage().getImplementationVersion();

    public static void main(String[] args) throws ParseException, IOException {
        LOGGER.info("Running Knowledgebase Generator v{}", VERSION);

        Options options = createOptions();
        CommandLine cmd = createCommandLine(args, options);

        if (!validInputForKnowledgebaseGeneration(cmd)) {
            printUsageAndExit(options);
        }

        //        String iClusionTrialTsv = cmd.getOptionValue(ICLUSION_TRIAL_TSV);
        //        List<IclusionTrial> trials = IclusionTrialFile.read(cmd.getOptionValue(ICLUSION_TRIAL_TSV));
        //        LOGGER.info("Read {} trials from {}", trials.size(), iClusionTrialTsv);
        //
        String viccJson = cmd.getOptionValue(VICC_JSON);
        List<ViccEntry> viccEntries = ViccJsonReader.readViccKnowledgebaseJsonFile(viccJson);
        LOGGER.info("Read {} VICC entries from {}", viccEntries.size(), viccJson);
        //
        //        String compassionateUseProgramsTsv = cmd.getOptionValue(COMPASSIONATE_USE_PROGRAMS_TSV);
        //        List<CompassionateUseProgram> compassionateUsePrograms = CompassionateUseProgramFile.read(compassionateUseProgramsTsv);
        //        LOGGER.info("Read {} compassionate use programs from {}", compassionateUsePrograms.size(), compassionateUseProgramsTsv);

        LOGGER.info("Convert VICC entries");

        String refFastaPath = "/data/common/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta";
        RefVersion refVersion = RefVersion.HG19;

        Transvar transvar = new Transvar(refFastaPath, refVersion);

        LOGGER.info("Generating known and actionable amps and dels");
        GeneratingCNV.generatingCNVs(viccEntries);

    }

    private static boolean validInputForKnowledgebaseGeneration(@NotNull CommandLine cmd) {
        return fileExists(cmd, ICLUSION_TRIAL_TSV) && fileExists(cmd, VICC_JSON) && fileExists(cmd, COMPASSIONATE_USE_PROGRAMS_TSV);
    }

    private static boolean fileExists(@NotNull CommandLine cmd, @NotNull String param) {
        if (!cmd.hasOption(param)) {
            LOGGER.warn("{} has to be provided", param);
            return false;
        } else if (!Files.exists(new File(cmd.getOptionValue(param)).toPath())) {
            LOGGER.warn("{} has to be an existing path", cmd.getOptionValue(param));
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
        options.addOption(COMPASSIONATE_USE_PROGRAMS_TSV, true, "compassionate use pgram input tsv");

        return options;
    }

    private static void printUsageAndExit(@NotNull Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Knowledgebase-Generator", options);
        System.exit(1);
    }
}
