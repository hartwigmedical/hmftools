package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.ReferenceData.TARGET_REGION_BED;

import java.io.File;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.purity.RunMode;
import com.hartwig.hmftools.common.variant.VariantTier;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class PurpleConfig
{
    public final String ReferenceId;
    public final String TumorId;

    public final String Version;
    public final String OutputDir;

    public final boolean RunDrivers;
    public final boolean DriversOnly;

    public final FittingConfig Fitting;
    public final SomaticFitConfig SomaticFitting;
    public final ChartConfig Charting;
    public final boolean TargetRegionsMode;
    public final Map<VariantTier,Integer> TierQualFilters;

    private boolean mIsValid;

    private static final String REF_SAMPLE = "reference";
    private static final String TUMOR_SAMPLE = "tumor";
    public static final String SAMPLE_DIR = "sample_dir";
    private static final String OUTPUT_DIRECTORY = "output_dir";
    private static final String AMBER = "amber";
    private static final String COBALT = "cobalt";

    public static String RUN_DRIVERS = "run_drivers";
    public static String DRIVERS_ONLY = "drivers_only";
    public static final String TIER_FILTERS = "tier_filters";

    public PurpleConfig(final String version, final CommandLine cmd)
    {
        mIsValid = true;

        Version = version;

        final StringJoiner missingJoiner = new StringJoiner(", ");

        String tumorId = cmd.getOptionValue(TUMOR_SAMPLE);
        ReferenceId = cmd.getOptionValue(REF_SAMPLE);

        TumorId = tumorId != null ? tumorId : ReferenceId;

        if(TumorId == null)
            missingJoiner.add(TUMOR_SAMPLE);

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
            PPL_LOGGER.error("unable to write directory " + OutputDir);
        }

        PPL_LOGGER.info("output directory: {}", OutputDir);

        Charting = new ChartConfig(cmd, OutputDir);
        Fitting = new FittingConfig(cmd);
        SomaticFitting = new SomaticFitConfig(cmd);
        TargetRegionsMode = cmd.hasOption(TARGET_REGION_BED);

        RunDrivers = cmd.hasOption(RUN_DRIVERS);
        DriversOnly = cmd.hasOption(DRIVERS_ONLY);

        PPL_LOGGER.info("reference({}) tumor({}) {}",
                ReferenceId != null ? ReferenceId : "NONE", TumorId != null ? TumorId : "NONE",
                TargetRegionsMode ? "running on target-regions only" : "");

        TierQualFilters = Maps.newHashMap();

        if(cmd.hasOption(TIER_FILTERS))
        {
            String[] tierFilterStrings = cmd.getOptionValue(TIER_FILTERS).split(";", -1);

            for(String tierFilter : tierFilterStrings)
            {
                String[] tierItems = tierFilter.split("=",-1);
                TierQualFilters.put(VariantTier.valueOf(tierItems[0]), Integer.parseInt(tierItems[1]));

                PPL_LOGGER.info("applying tier({}) qual({}) filter", tierItems[0], tierItems[1]);
            }
        }
    }

    public boolean isValid() { return mIsValid; }

    public boolean tumorOnlyMode() { return ReferenceId == null; }
    public boolean germlineMode() { return TumorId.equals(ReferenceId); }

    public boolean runTumor() { return !germlineMode(); }
    public boolean runGermline() { return !tumorOnlyMode(); }

    public boolean fitWithSomatics() { return !tumorOnlyMode() && !germlineMode() && !TargetRegionsMode; }

    public RunMode runMode()
    {
        return tumorOnlyMode() ? RunMode.TUMOR : (germlineMode() ? RunMode.GERMLINE : RunMode.TUMOR_GERMLINE);
    }

    public static void addOptions(@NotNull Options options)
    {
        options.addOption(
                REF_SAMPLE, true,
                "Name of the reference sample. This should correspond to the value used in Amber and Cobalt");

        options.addOption(
                TUMOR_SAMPLE, true,
                "Name of the tumor sample. This should correspond to the value used in Amber and Cobalt");

        options.addOption(
                OUTPUT_DIRECTORY, true,
                "Path to the output directory. Required if <run_dir> not set, otherwise defaults to run_dir/purple/");

        options.addOption(
                SAMPLE_DIR, true,
                "Path to the sample's directory where expect to find Cobalt, Amber, Gripss etc directories");

        options.addOption(
                COBALT, true,
                "Path to Cobalt output directory. Required if <run_dir> not set, otherwise defaults to <run_dir>/cobalt");

        options.addOption(
                AMBER, true,
                "Path to Amber output directory. Required if <run_dir> not set, otherwise defaults to <run_dir>/amber");

        options.addOption(RUN_DRIVERS, false, "Run driver routine");
        options.addOption(DRIVERS_ONLY, false, "Only run the driver routine");
        options.addOption(TIER_FILTERS, true, "Variant qual filters by tier, format: TIER_A=QUAL;TIER_A=QUAL etc");

        addDatabaseCmdLineArgs(options);
        FittingConfig.addOptions(options);
        SomaticFitConfig.addOptions(options);
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
