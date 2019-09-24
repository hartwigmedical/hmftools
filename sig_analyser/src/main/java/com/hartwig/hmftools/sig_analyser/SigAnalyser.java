package com.hartwig.hmftools.sig_analyser;

import java.sql.SQLException;

import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.sig_analyser.buckets.BucketAnalyser;
import com.hartwig.hmftools.sig_analyser.common.CosineSim;
import com.hartwig.hmftools.sig_analyser.nmf.NmfConfig;
import com.hartwig.hmftools.sig_analyser.nmf.NmfManager;
import com.hartwig.hmftools.sig_analyser.sim.SampleSimulator;
import com.hartwig.hmftools.sig_analyser.sim.SimConfig;
import com.hartwig.hmftools.common.utils.GenericDataLoader;
import com.hartwig.hmftools.common.utils.GenericDataCollection;
import com.hartwig.hmftools.sig_analyser.loaders.SigSnvLoader;

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

public class SigAnalyser
{
    private static final String RUN_CSS = "run_cosine_sim";
    private static final String RUN_NMF = "run_nmf";
    private static final String RUN_BA = "run_buckets";
    private static final String RUN_SAMPLE_SIM = "run_sim";
    private static final String LOAD_SNVS = "load_snvs";

    private static final String GENERIC_INPUT_FILE = "gen_input_file";
    public static final String SAMPLE_IDS = "sample_ids";
    private static final String LOG_DEBUG = "log_debug";

    public static final String OUTPUT_DIR = "output_dir";
    public static final String OUTPUT_FILE_ID = "output_file_id";

    private static final Logger LOGGER = LogManager.getLogger(SigAnalyser.class);

    public static void main(@NotNull final String[] args) throws ParseException {

        Options options = createBasicOptions();

        // allow components to add their own arg lists
        SimConfig.addCmdLineArgs(options);
        NmfConfig.addCmdLineArgs(options);
        BucketAnalyser.addCmdLineArgs(options);
        SigSnvLoader.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        if(cmd.hasOption(LOAD_SNVS))
        {
            try
            {
                final DatabaseAccess dbAccess = cmd.hasOption(DB_URL) ? databaseAccess(cmd) : null;

                SigSnvLoader snvLoader = new SigSnvLoader();

                if(snvLoader.initialise(dbAccess, cmd))
                {
                    snvLoader.loadData();
                }

                LOGGER.info("variant download complete");
            }
            catch(SQLException e)
            {
                LOGGER.error("DB connection failed: {}", e.toString());
            }

            return;
        }

        LOGGER.info("starting signature analyser");

        GenericDataCollection collection = GenericDataLoader.loadFile(cmd.getOptionValue(GENERIC_INPUT_FILE));

        if(cmd.hasOption(RUN_CSS))
        {

            CosineSim cosineSim = new CosineSim(cmd.getOptionValue(OUTPUT_DIR));
            cosineSim.calcCosineSimilarities(collection.getFieldNames(), collection.getData(), 0.8);
        }

        if(cmd.hasOption(RUN_NMF))
        {

            NmfManager nmfManager = new NmfManager();
            nmfManager.initialise(collection, cmd);
            nmfManager.run();
            // nmfManager.runTests();
        }

        if(cmd.hasOption(RUN_BA))
        {
            BucketAnalyser bucketAnalyser = new BucketAnalyser();
            if(!bucketAnalyser.initialise(collection, cmd))
            {
                LOGGER.info("init failed");
                return;
            }

            bucketAnalyser.run();
        }

        if(cmd.hasOption(RUN_SAMPLE_SIM))
        {
            SampleSimulator sampleSimulator = new SampleSimulator();
            sampleSimulator.initialise(cmd);
            // sampleSimulator.runTests();
            sampleSimulator.run();
        }

        LOGGER.info("analysis complete");
    }

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    @NotNull
    private static Options createBasicOptions()
    {
        Options options = new Options();
        options.addOption(GENERIC_INPUT_FILE, true, "Path to the main input file.");
        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(OUTPUT_FILE_ID, true, "Output file ID");

        options.addOption(RUN_CSS, false, "Run cosine similiarities");
        options.addOption(RUN_NMF, false, "Run NMF");
        options.addOption(RUN_BA, false, "Run bucket analysis");
        options.addOption(RUN_SAMPLE_SIM, false, "Generate simulated sample counts");
        options.addOption(LOAD_SNVS, false, "Create sample bucket counts for SNVs from DB");
        options.addOption(SAMPLE_IDS, true, "Optional: restrict to sample list");

        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull final CommandLine cmd) throws SQLException
    {
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }


}