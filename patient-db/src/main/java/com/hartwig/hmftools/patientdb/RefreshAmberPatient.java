package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.sql.SQLException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.amber.AmberMapping;
import com.hartwig.hmftools.common.amber.AmberMappingFactory;
import com.hartwig.hmftools.common.amber.AmberPatient;
import com.hartwig.hmftools.common.amber.AmberPatientFactory;
import com.hartwig.hmftools.common.amber.AmberSample;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class RefreshAmberPatient {

    private static final Logger LOGGER = LogManager.getLogger(RefreshAmberPatient.class);

    public static void main(@NotNull String[] args) throws ParseException, SQLException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        try (DatabaseAccess dbAccess = databaseAccess(cmd)) {
            LOGGER.info("Reading sample data");
            List<AmberPatient> previousPatients = dbAccess.readAmberPatients();
            List<AmberSample> allSamples = dbAccess.readAmberSamples();
            List<AmberMapping> allMappings = Lists.newArrayList();

            for (int i = 0; i < allSamples.size(); i++) {
                AmberSample victim = allSamples.get(i);

                LOGGER.info("Processing " + (i + 1) + " of " + allSamples.size() + ": " + victim.sampleId());
                for (int j = i + 1; j < allSamples.size(); j++) {
                    AmberMapping match = AmberMappingFactory.create(victim, allSamples.get(j));
                    if (match.likelihood() > 0.8) {
                        allMappings.add(match);
                    }
                }
            }

            AmberPatientFactory patientFactory = new AmberPatientFactory(previousPatients, allMappings);
            List<AmberPatient> patients = allSamples.stream().map(patientFactory::createPatient).collect(Collectors.toList());

            LOGGER.info("Truncating AMBER mappings and patients");
            dbAccess.truncateAmberPatients();
            dbAccess.truncateAmberMappings();

            LOGGER.info("Writing mapping data");
            dbAccess.writeAmberMapping("dummy", allMappings);
            dbAccess.writeAmberPatients(patients);
        }

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        addDatabaseCmdLineArgs(options);
        return options;
    }
}
