package com.hartwig.hmftools.patientdb.amber;

import static com.hartwig.hmftools.patientdb.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class RefreshAmberPatient
{
    public static void main(@NotNull String[] args) throws ParseException, SQLException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        addDatabaseCmdLineArgs(configBuilder, true);
        configBuilder.checkAndParseCommandLine(args);

        try(DatabaseAccess dbAccess = databaseAccess(configBuilder))
        {
            LOGGER.info("Reading sample data");
            List<AmberPatient> previousPatients = dbAccess.readAmberPatients();
            List<AmberSample> allSamples = dbAccess.readAmberSamples();
            List<AmberMapping> allMappings = new ArrayList<>();

            for(int i = 0; i < allSamples.size(); i++)
            {
                AmberSample victim = allSamples.get(i);

                LOGGER.info("Processing " + (i + 1) + " of " + allSamples.size() + ": " + victim.sampleId());
                for(int j = i + 1; j < allSamples.size(); j++)
                {
                    AmberMapping match = AmberMappingFactory.create(victim, allSamples.get(j));
                    if(match.likelihood() > 0.8)
                    {
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
}
