package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

import java.io.File;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class PurpleConfig
{
    public final String ReferenceId;
    public final String TumorId;

    public final String Version;
    public final String OutputDir;

    public final boolean TumorOnlyMode;
    public final boolean RunDrivers;
    public final boolean DriversOnly;

    public final FittingConfig Fitting;
    public final SomaticFitConfig SomaticFitting;
    public final MiscConfig Misc;
    public final ChartConfig Charting;

    private boolean mIsValid;

    private static final String REF_SAMPLE = "reference";
    private static final String TUMOR_SAMPLE = "tumor";
    public static final String SAMPLE_DIR = "sample_dir";
    private static final String OUTPUT_DIRECTORY = "output_dir";
    private static final String AMBER = "amber";
    private static final String COBALT = "cobalt";

    public static String RUN_DRIVERS = "run_drivers";
    public static String DRIVERS_ONLY = "drivers_only";

    public PurpleConfig(final String version, final CommandLine cmd)
    {
        mIsValid = true;

        Version = version;

        final StringJoiner missingJoiner = new StringJoiner(", ");

        RunDrivers = cmd.hasOption(RUN_DRIVERS);
        DriversOnly = cmd.hasOption(DRIVERS_ONLY);

        TumorId = parameter(cmd, TUMOR_SAMPLE, missingJoiner);
        ReferenceId = cmd.getOptionValue(REF_SAMPLE);
        TumorOnlyMode = ReferenceId == null;

        String sampleDir = "";

        if(cmd.hasOption(SAMPLE_DIR))
        {
            sampleDir = checkAddDirSeparator(cmd.getOptionValue(SAMPLE_DIR));
            OutputDir = checkAddDirSeparator(sampleDir + cmd.getOptionValue(OUTPUT_DIRECTORY, ""));
        }
        else
        {
            OutputDir = checkAddDirSeparator(parameter(cmd, OUTPUT_DIRECTORY, missingJoiner));
        }

        parameter(cmd, ENSEMBL_DATA_DIR, missingJoiner);

        final String missing = missingJoiner.toString();

        if(!missing.isEmpty())
        {
            mIsValid = false;
            PPL_LOGGER.error("Missing the following parameters: " + missing);
        }

        final File outputDir = new File(OutputDir);
        if(!outputDir.exists() && !outputDir.mkdirs())
        {
            mIsValid = false;
            PPL_LOGGER.error("Unable to write directory " + OutputDir);
        }

        if(TumorOnlyMode)
        {
            PPL_LOGGER.info("Tumor Sample: {}", TumorId);
        }
        else
        {
            PPL_LOGGER.info("Reference Sample: {}, Tumor Sample: {}", ReferenceId, TumorId);
        }

        PPL_LOGGER.info("Output Directory: {}", OutputDir);

        Charting = new ChartConfig(cmd, OutputDir);
        Fitting = new FittingConfig(cmd);
        SomaticFitting = new SomaticFitConfig(cmd);
        Misc = new MiscConfig(cmd);
    }

    public boolean isValid() { return mIsValid; }

    public static void addOptions(@NotNull Options options)
    {
        options.addOption(REF_SAMPLE, true, "Name of the reference sample. This should correspond to the value used in AMBER and COBALT.");
        options.addOption(TUMOR_SAMPLE, true, "Name of the tumor sample. This should correspond to the value used in AMBER and COBALT.");

        options.addOption(OUTPUT_DIRECTORY,
                true,
                "Path to the output directory. Required if <run_dir> not set, otherwise defaults to run_dir/purple/");

        options.addOption(SAMPLE_DIR,
                true,
                "Path to the sample's directory where expect to find cobalt, amber, gridss etc directories");

        options.addOption(COBALT,
                true,
                "Path to COBALT output directory. Required if <run_dir> not set, otherwise defaults to <run_dir>/cobalt.");

        options.addOption(AMBER,
                true,
                "Path to AMBER output directory. Required if <run_dir> not set, otherwise defaults to <run_dir>/amber");

        options.addOption(RUN_DRIVERS, false, "Run driver routine");
        options.addOption(DRIVERS_ONLY, false, "Only run the driver routine");

        addDatabaseCmdLineArgs(options);
        FittingConfig.addOptions(options);
        SomaticFitConfig.addOptions(options);
        MiscConfig.addOptions(options);
        ReferenceData.addOptions(options);
        ChartConfig.addOptions(options);
        SampleDataFiles.addOptions(options);
    }

    private static String parameter(final CommandLine cmd, final String parameter, final StringJoiner missing)
    {
        final String value = cmd.getOptionValue(parameter);
        if(value == null)
        {
            missing.add(parameter);
            return "";
        }
        return value;
    }
}
