package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_CHROMOSOMES;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificChromsomesOrRegions;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.purple.RunMode;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.VariantTier;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;

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
    public final List<String> SpecificChromosomes;
    public final List<ChrBaseRegion> SpecificRegions;

    private boolean mIsValid;

    private static final String REF_SAMPLE = "reference";
    private static final String TUMOR_SAMPLE = "tumor";
    public static final String SAMPLE_DIR = "sample_dir";
    private static final String OUTPUT_DIRECTORY = "output_dir";
    private static final String AMBER = "amber";
    private static final String COBALT = "cobalt";

    public static String DRIVERS_ONLY = "drivers_only";
    public static String FILTER_SOMATICS_ON_GENE = "filter_somatics_on_gene";
    public static final String TIER_FILTERS = "tier_filters";
    public static final String WRITE_ALL_SOMATICS = "write_all_somatics";

    public PurpleConfig(final String version, final ConfigBuilder configBuilder)
    {
        mIsValid = true;

        Version = version;

        String tumorId = configBuilder.getValue(TUMOR_SAMPLE);
        ReferenceId = configBuilder.getValue(REF_SAMPLE);
        TumorId = tumorId != null ? tumorId : ReferenceId;

        if(TumorId == null && ReferenceId == null)
        {
            mIsValid = false;
        }

        String outputDir = configBuilder.getValue(OUTPUT_DIRECTORY);
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

        final File outputPath = new File(OutputDir);
        if(!outputPath.exists() && !outputPath.mkdirs())
        {
            mIsValid = false;
            PPL_LOGGER.error("unable to write directory " + OutputDir);
        }

        PPL_LOGGER.info("output directory: {}", OutputDir);

        SampleFiles = new SampleDataFiles(configBuilder, TumorId);

        Charting = new ChartConfig(configBuilder, OutputDir);
        Fitting = new FittingConfig(configBuilder);
        SomaticFitting = new SomaticFitConfig(configBuilder);
        TargetRegionsMode = configBuilder.hasValue(TARGET_REGIONS_BED);
        Threads = parseThreads(configBuilder);

        RunDrivers = DriverGenePanelConfig.isConfigured(configBuilder);
        DriversOnly = configBuilder.hasFlag(DRIVERS_ONLY);
        FilterSomaticsOnGene = configBuilder.hasFlag(FILTER_SOMATICS_ON_GENE);
        WriteAllSomatics = configBuilder.hasFlag(WRITE_ALL_SOMATICS);

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

        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();

        try
        {
            loadSpecificChromsomesOrRegions(configBuilder, SpecificChromosomes, SpecificRegions, PPL_LOGGER);
            Collections.sort(SpecificRegions);
        }
        catch(ParseException e)
        {
            PPL_LOGGER.error("invalid specific regions({}) chromosomes({}) config",
                    configBuilder.getValue(SPECIFIC_REGIONS), configBuilder.getValue(SPECIFIC_CHROMOSOMES));
            mIsValid = false;
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

    public static void addOptions(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(
                REF_SAMPLE, "Name of the reference sample. This should correspond to the value used in Amber and Cobalt");

        configBuilder.addConfigItem(
                TUMOR_SAMPLE,
                "Name of the tumor sample. This should correspond to the value used in Amber and Cobalt");

        configBuilder.addConfigItem(
                OUTPUT_DIRECTORY, true,
                "Path to the output directory. If <sample_dir> is set, then is sample_dir/output_dir/",
                "");

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
        if(!SpecificRegions.isEmpty())
            return SpecificRegions.stream().noneMatch(x -> x.containsPosition(chromosome, position));

        if(!SpecificChromosomes.isEmpty())
            return !SpecificChromosomes.contains(chromosome);

        return false;
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
