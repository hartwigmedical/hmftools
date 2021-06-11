package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.AnnotatedVirusFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadVirusInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(LoadVirusInterpreter.class);

    private static final String SAMPLE = "sample";
    private static final String VIRUS_ANNOTATION_TSV = "virus_annotation_tsv";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String sample = cmd.getOptionValue(SAMPLE);
        String virusAnnotationTsv = cmd.getOptionValue(VIRUS_ANNOTATION_TSV);

        if (Utils.anyNull(sample, virusAnnotationTsv)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-DB - Load VirusInterpreter Data", options);
            System.exit(1);
        }

        DatabaseAccess dbWriter = databaseAccess(cmd);

        LOGGER.info("Reading virus annotation TSV {}", virusAnnotationTsv);
        List<AnnotatedVirus> virusAnnotations = AnnotatedVirusFile.read(virusAnnotationTsv);
        LOGGER.info(" Read {} virus annotation", virusAnnotations.size());

        LOGGER.info("Writing virus annotations into database for {}", sample);
        dbWriter.writeVirusInterpreter(sample, virusAnnotations);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(SAMPLE, true, "Sample for which we are going to load the virus breakends");
        options.addOption(VIRUS_ANNOTATION_TSV, true, "Path towards the virus annotations TSV file");

        addDatabaseCmdLineArgs(options);

        return options;
    }
}
