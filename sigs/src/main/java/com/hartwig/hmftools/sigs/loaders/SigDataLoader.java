package com.hartwig.hmftools.sigs.loaders;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;
import static com.hartwig.hmftools.sigs.common.CommonUtils.formOutputFilename;

import java.sql.SQLException;

import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class SigDataLoader
{
    // config
    private static final String LOAD_SNVS = "load_snvs";
    private static final String LOAD_MNVS = "load_mnvs";
    private static final String LOAD_INDELS = "load_indels";

    private static final Logger LOGGER = LogManager.getLogger(SigDataLoader.class);

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = createBasicOptions();

        DataLoaderConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        final DataLoaderConfig config = new DataLoaderConfig(cmd);

        try
        {
            final DatabaseAccess dbAccess = databaseAccess(cmd);

            config.loadSampleIds(dbAccess);

            if(cmd.hasOption(LOAD_SNVS))
            {
                SigSnvLoader snvLoader = new SigSnvLoader(config.Filters);
                snvLoader.setSampleIds(config.SampleIds);

                if(!config.PositionBucketSizes.isEmpty())
                    snvLoader.initialisePositionFrequencies(config.OutputDir, config.PositionBucketSizes);

                snvLoader.loadData(dbAccess, null, true);

                final String fileId = config.SampleIds.size() == 1 ? config.SampleIds.get(0) + "." + "sample_counts" : "sample_counts";
                final String filename = formOutputFilename(config.OutputDir, config.OutputFileId, fileId);

                snvLoader.writeSampleCounts(filename);
            }

            if(cmd.hasOption(LOAD_MNVS))
            {
                SigMnvLoader snvLoader = new SigMnvLoader(config);
                snvLoader.loadData(dbAccess);
            }

            if(cmd.hasOption(LOAD_INDELS))
            {
                SigIndelLoader snvLoader = new SigIndelLoader(config);
                snvLoader.loadData(dbAccess);
            }
        }
        catch(SQLException e)
        {
            LOGGER.error("DB connection failed: {}", e.toString());
        }

        LOGGER.info("data load complete");
    }

    @NotNull
    private static Options createBasicOptions()
    {
        Options options = new Options();

        options.addOption(LOAD_SNVS, false, "Create sample bucket counts for SNVs from DB");
        options.addOption(LOAD_MNVS, false, "Create sample bucket counts for MNVs from DB");
        options.addOption(LOAD_INDELS, false, "Create sample bucket counts for INDELs from DB");

        addDatabaseCmdLineArgs(options);

        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
