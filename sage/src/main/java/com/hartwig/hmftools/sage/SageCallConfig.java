package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAMS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_IDS_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.sage.SageCommon.SAMPLE_DELIM;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConfig.registerCommonConfig;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.util.Strings;

public class SageCallConfig
{
    public final SageConfig Common;

    public final List<String> TumorIds;
    public final List<String> TumorBams;
    public final String HighConfidenceBed;
    public final String CoverageBed;
    public final String PanelBed;
    public final String Hotspots;
    public final boolean PanelOnly;

    private final String mResourceDir;

    private static final String COVERAGE_BED = "coverage_bed";
    private static final String RESOURCE_DIR = "resource_dir";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";
    private static final String PANEL_BED = "panel_bed";
    private static final String HOTSPOTS = "hotspots";
    private static final String PANEL_ONLY = "panel_only";

    public SageCallConfig(final String version, final ConfigBuilder configBuilder)
    {
        Common = new SageConfig(version, configBuilder);

        TumorIds = Lists.newArrayList();
        if(configBuilder.hasValue(TUMOR))
        {
            TumorIds.addAll(Arrays.asList(configBuilder.getValue(TUMOR).split(SAMPLE_DELIM)));
        }

        TumorBams = Lists.newArrayList();

        if(configBuilder.hasValue(TUMOR_BAM))
        {
            Arrays.stream(configBuilder.getValue(TUMOR_BAM, Strings.EMPTY).split(SAMPLE_DELIM))
                    .forEach(x -> TumorBams.add(Common.SampleDataDir + x));
        }

        mResourceDir = checkAddDirSeparator(configBuilder.getValue(RESOURCE_DIR, ""));
        PanelBed = getReferenceFile(configBuilder, PANEL_BED);
        CoverageBed = getReferenceFile(configBuilder, COVERAGE_BED);
        HighConfidenceBed = getReferenceFile(configBuilder, HIGH_CONFIDENCE_BED);
        Hotspots = getReferenceFile(configBuilder, HOTSPOTS);

        PanelOnly = configBuilder.hasFlag(PANEL_ONLY);
    }

    public boolean isValid()
    {
        if(!Common.isValid())
            return false;

        if(TumorIds.size() != TumorBams.size())
        {
            SG_LOGGER.error("Each tumor sample must have matching bam");
            return false;
        }

        for(String tumorBam : TumorBams)
        {
            if(!new File(tumorBam).exists())
            {
                SG_LOGGER.error("Unable to locate tumor bam({})", tumorBam);
                return false;
            }
        }

        if(TumorIds.isEmpty())
        {
            SG_LOGGER.error("At least one tumor must be supplied");
            return false;
        }

        return true;
    }

    private String getReferenceFile(final ConfigBuilder configBuilder, final String config)
    {
        if(!configBuilder.hasValue(config))
            return "";

        if(mResourceDir.isEmpty())
            return configBuilder.getValue(config);

        return mResourceDir + configBuilder.getValue(config);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(TUMOR, true, TUMOR_IDS_DESC);
        configBuilder.addConfigItem(TUMOR_BAM, true, TUMOR_BAMS_DESC);

        configBuilder.addPath(RESOURCE_DIR, false, "Path to Sage resource files");
        configBuilder.addPrefixedPath(HIGH_CONFIDENCE_BED, false, "High confidence regions bed file", RESOURCE_DIR);
        configBuilder.addPrefixedPath(PANEL_BED, false, "Panel regions bed file", RESOURCE_DIR);
        configBuilder.addPrefixedPath(HOTSPOTS, false, "Hotspots", RESOURCE_DIR);
        configBuilder.addPrefixedPath(COVERAGE_BED, false, "Coverage is calculated for optionally supplied bed", RESOURCE_DIR);
        configBuilder.addFlag(PANEL_ONLY, "Only examine panel for variants");

        registerCommonConfig(configBuilder);
        addEnsemblDir(configBuilder);
    }

    public SageCallConfig()
    {
        Common = new SageConfig(false);
        TumorIds = Lists.newArrayList();
        TumorBams = Lists.newArrayList();
        HighConfidenceBed = "highConf";
        CoverageBed = "coverage";
        PanelBed = "panel";
        Hotspots = "hotspots";
        PanelOnly = false;
        mResourceDir = "";
    }
}
