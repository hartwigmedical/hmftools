package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.sql.SQLException;

import com.hartwig.hmftools.common.genome.region.CanonicalTranscriptFactory;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadCanonicalTranscripts {

    private static final Logger LOGGER = LogManager.getLogger(LoadCanonicalTranscripts.class);

    public static void main(@NotNull String[] args) throws ParseException, SQLException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);
        DatabaseAccess dbAccess = databaseAccess(cmd);

        LOGGER.info("Persisting transcripts to database");
        dbAccess.writeCanonicalTranscripts("GRCh37", CanonicalTranscriptFactory.create37());
        dbAccess.writeCanonicalTranscripts("GRCh38", CanonicalTranscriptFactory.create38());

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        addDatabaseCmdLineArgs(options);
        return options;
    }
}
