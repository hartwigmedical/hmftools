package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.File;

import com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.purple.RunMode;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class PurpleConfig
{
    public final String ReferenceId;
    public final String TumorId;

    public final String Version;
    public final String OutputDir;

    public final boolean RunDrivers;

    public final SampleDataFiles SampleFiles;
    public final FittingConfig Fitting;
    public final SomaticFitConfig SomaticFitting;
    public final ChartConfig Charting;
    public final boolean TargetRegionsMode;
    public final boolean IgnorePlotErrors;
    public final int Threads;

    // debug only
    public final boolean FilterSomaticsOnGene;
    public final boolean WriteAllSomatics;
    public final SpecificRegions SpecificChrRegions;

    private boolean mIsValid;

    public static final String SAMPLE_DIR = "sample_dir";

    public static String FILTER_SOMATICS_ON_GENE = "filter_somatics_on_gene";
    public static final String WRITE_ALL_SOMATICS = "write_all_somatics";
    public static final String IGNORE_PLOT_ERRORS = "ignore_plot_errors";

    public PurpleConfig(final String version, final ConfigBuilder configBuilder)
    {
        mIsValid = true;

        Version = version;

        String tumorId = configBuilder.getValue(TUMOR);
        ReferenceId = configBuilder.getValue(REFERENCE);
        TumorId = tumorId != null ? tumorId : ReferenceId;

        if(TumorId == null)
        {
            mIsValid = false;
        }

        String outputDir = configBuilder.getValue(OUTPUT_DIR);
        String sampleDir;

        if(configBuilder.hasValue(SAMPLE_DIR))
        {
            sampleDir = checkAddDirSeparator(configBuilder.getValue(SAMPLE_DIR));
            OutputDir = checkAddDirSeparator(sampleDir + outputDir);
        }
        else
        {
            OutputDir = checkAddDirSeparator(outputDir);
        }

        SampleFiles = new SampleDataFiles(configBuilder, TumorId);

        Charting = new ChartConfig(configBuilder, OutputDir);

        if(OutputDir != null)
        {
            mIsValid &= createDirectory(OutputDir);

            if(!Charting.Disabled)
            {
                if(Charting.CircosBinary != null)
                {
                    mIsValid &= createDirectory(Charting.CircosDirectory);
                }

                mIsValid &= createDirectory(Charting.PlotDirectory);
            }
        }

        TargetRegionsMode = configBuilder.hasValue(TARGET_REGIONS_BED);
        Fitting = new FittingConfig(configBuilder, TargetRegionsMode);
        SomaticFitting = new SomaticFitConfig(configBuilder);
        Threads = parseThreads(configBuilder);

        RunDrivers = DriverGenePanelConfig.isConfigured(configBuilder);
        FilterSomaticsOnGene = configBuilder.hasFlag(FILTER_SOMATICS_ON_GENE);
        IgnorePlotErrors = configBuilder.hasFlag(IGNORE_PLOT_ERRORS);
        WriteAllSomatics = configBuilder.hasFlag(WRITE_ALL_SOMATICS);

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        if(SpecificChrRegions == null)
        {
            mIsValid = false;
        }
    }

    public boolean isValid()
    {
        return mIsValid;
    }
    public boolean tumorOnlyMode()
    {
        return ReferenceId == null;
    }
    public boolean germlineMode()
    {
        return TumorId.equals(ReferenceId);
    }
    public boolean runTumor()
    {
        return !germlineMode();
    }
    public boolean runGermline()
    {
        return !tumorOnlyMode();
    }
    public RunMode runMode()
    {
        return tumorOnlyMode() ? RunMode.TUMOR : (germlineMode() ? RunMode.GERMLINE : RunMode.TUMOR_GERMLINE);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(REFERENCE, REFERENCE_DESC);
        configBuilder.addConfigItem(TUMOR, TUMOR_DESC);

        configBuilder.addConfigItem(
                OUTPUT_DIR, false,
                "Path to the output directory. If <sample_dir> is set, then is sample_dir/output_dir/.");

        configBuilder.addFlag(WRITE_ALL_SOMATICS, "Write all variants regardless of filters");
        configBuilder.addFlag(FILTER_SOMATICS_ON_GENE, "Only load and enrich somatic variants with a gene impact");
        configBuilder.addFlag(IGNORE_PLOT_ERRORS, "Run to completion if plotting fails");

        FittingConfig.addConfig(configBuilder);
        SomaticFitConfig.addConfig(configBuilder);
        ReferenceData.addConfig(configBuilder);
        ChartConfig.addConfig(configBuilder);
        SampleDataFiles.addConfig(configBuilder);
        addThreadOptions(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);
    }

    public boolean excludeOnSpecificRegion(final String chromosome, final int position)
    {
        if(!SpecificChrRegions.hasFilters())
        {
            return false;
        }

        if(!SpecificChrRegions.Regions.isEmpty())
        {
            return SpecificChrRegions.excludePosition(chromosome, position);
        }

        return SpecificChrRegions.excludeChromosome(chromosome);
    }

    private boolean createDirectory(final String dir)
    {
        final File output = new File(dir);
        if(!output.exists() && !output.mkdirs())
        {
            PPL_LOGGER.error("unable to create chart directory " + dir);
            return false;
        }

        return true;
    }
}
