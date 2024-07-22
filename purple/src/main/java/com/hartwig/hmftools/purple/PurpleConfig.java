package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PURPLE_DIR;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.File;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.purple.RunMode;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VariantTier;

public class PurpleConfig
{
    public final String ReferenceId;
    public final String TumorId;

    public final String Version;
    public final String OutputDir;

    public final boolean RunDrivers;
    public final boolean DriversOnly;

    public final SampleDataFiles SampleFiles;
    public final FittingConfig Fitting;
    public final SomaticFitConfig SomaticFitting;
    public final ChartConfig Charting;
    public final boolean TargetRegionsMode;
    public final Map<VariantTier,Integer> TierQualFilters;
    public final int Threads;

    // debug only
    public final boolean FilterSomaticsOnGene;
    public final boolean WriteAllSomatics;
    public final boolean UseGridssSVs;
    public final SpecificRegions SpecificChrRegions;

    private boolean mIsValid;

    public static final String SAMPLE_DIR = "sample_dir";

    public static String DRIVERS_ONLY = "drivers_only";
    public static String FILTER_SOMATICS_ON_GENE = "filter_somatics_on_gene";
    public static final String TIER_FILTERS = "tier_filters";
    public static final String WRITE_ALL_SOMATICS = "write_all_somatics";

    public PurpleConfig(final String version, final ConfigBuilder configBuilder)
    {
        mIsValid = true;

        Version = version;

        String tumorId = configBuilder.getValue(TUMOR);
        ReferenceId = configBuilder.getValue(REFERENCE);
        TumorId = tumorId != null ? tumorId : ReferenceId;

        if(TumorId == null && ReferenceId == null)
        {
            mIsValid = false;
        }

        String outputDir = configBuilder.getValue(OUTPUT_DIR);
        String sampleDir = "";

        if(configBuilder.hasValue(SAMPLE_DIR))
        {
            sampleDir = checkAddDirSeparator(configBuilder.getValue(SAMPLE_DIR));
            OutputDir = checkAddDirSeparator(sampleDir + outputDir);
        }
        else
        {
            OutputDir = checkAddDirSeparator(outputDir);
        }

        mIsValid &= createDirectory(OutputDir);

        PPL_LOGGER.info("output directory: {}", OutputDir);

        SampleFiles = new SampleDataFiles(configBuilder, TumorId);

        Charting = new ChartConfig(configBuilder, OutputDir);

        if(!Charting.Disabled)
        {
            if(Charting.CircosBinary != null)
                mIsValid &= createDirectory(Charting.CircosDirectory);

            mIsValid &= createDirectory(Charting.PlotDirectory);
        }

        TargetRegionsMode = configBuilder.hasValue(TARGET_REGIONS_BED);
        Fitting = new FittingConfig(configBuilder, TargetRegionsMode);
        SomaticFitting = new SomaticFitConfig(configBuilder);
        Threads = parseThreads(configBuilder);

        RunDrivers = DriverGenePanelConfig.isConfigured(configBuilder);
        DriversOnly = configBuilder.hasFlag(DRIVERS_ONLY);
        FilterSomaticsOnGene = configBuilder.hasFlag(FILTER_SOMATICS_ON_GENE);
        WriteAllSomatics = configBuilder.hasFlag(WRITE_ALL_SOMATICS);
        UseGridssSVs = SampleFiles.usesGridssSVs();

        if(UseGridssSVs)
        {
            PPL_LOGGER.info("using deprecated Gridss/Gripss VCFs");
        }

        PPL_LOGGER.info("reference({}) tumor({}) {}",
                ReferenceId != null ? ReferenceId : "NONE", TumorId != null ? TumorId : "NONE",
                TargetRegionsMode ? "running on target-regions only" : "");

        TierQualFilters = Maps.newHashMap();

        if(configBuilder.hasFlag(TIER_FILTERS))
        {
            String[] tierFilterStrings = configBuilder.getValue(TIER_FILTERS).split(";", -1);

            for(String tierFilter : tierFilterStrings)
            {
                String[] tierItems = tierFilter.split("=",-1);
                TierQualFilters.put(VariantTier.valueOf(tierItems[0]), Integer.parseInt(tierItems[1]));

                PPL_LOGGER.info("applying tier({}) qual({}) filter", tierItems[0], tierItems[1]);
            }
        }

        SpecificChrRegions = SpecificRegions.from(configBuilder);

        if(SpecificChrRegions == null)
            mIsValid = false;
    }

    public boolean isValid() { return mIsValid; }

    public boolean tumorOnlyMode() { return ReferenceId == null; }
    public boolean germlineMode() { return TumorId.equals(ReferenceId); }

    public boolean runTumor() { return !germlineMode(); }
    public boolean runGermline() { return !tumorOnlyMode(); }

    public boolean fitWithSomatics() { return !germlineMode(); }

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
                "Path to the output directory. If <sample_dir> is set, then is sample_dir/output_dir/. Default 'purple'",
                PURPLE_DIR);

        configBuilder.addFlag(DRIVERS_ONLY, "Only run the driver routine");
        configBuilder.addFlag(WRITE_ALL_SOMATICS, "Write all variants regardless of filters");
        configBuilder.addFlag(FILTER_SOMATICS_ON_GENE, "Only load and enrich somatic variants with a gene impact");
        configBuilder.addConfigItem(TIER_FILTERS, "Variant qual filters by tier, format: TIER_A=QUAL;TIER_A=QUAL etc");

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
            return false;

        if(!SpecificChrRegions.Regions.isEmpty())
            return SpecificChrRegions.excludePosition(chromosome, position);

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
