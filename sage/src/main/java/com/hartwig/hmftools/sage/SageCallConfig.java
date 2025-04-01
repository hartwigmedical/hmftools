package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAMS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_IDS_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.sage.SageCommon.SAMPLE_DELIM;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.sage.tinc.TincConfig;

import org.apache.logging.log4j.util.Strings;

public class SageCallConfig
{
    public final SageConfig Common;

    public final List<String> TumorIds;
    public final List<String> TumorBams;
    public final String HighConfidenceBed;
    public final String PanelBed;
    public final String Hotspots;
    public final boolean PanelOnly;
    public final boolean RunTinc;

    private final String mResourceDir;

    private static final String RESOURCE_DIR = "resource_dir";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";
    private static final String PANEL_BED = "panel_bed";
    private static final String HOTSPOTS = "hotspots";
    private static final String PANEL_ONLY = "panel_only";

    public static final String RUN_TINC = "run_tinc";

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
        HighConfidenceBed = getReferenceFile(configBuilder, HIGH_CONFIDENCE_BED);
        Hotspots = getReferenceFile(configBuilder, HOTSPOTS);

        PanelOnly = configBuilder.hasFlag(PANEL_ONLY);
        RunTinc = configBuilder.hasFlag(RUN_TINC);
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
        SageConfig.registerCommonConfig(configBuilder);

        configBuilder.addPath(RESOURCE_DIR, false, "Path to Sage resource files");
        configBuilder.addPrefixedPath(HIGH_CONFIDENCE_BED, false, "High confidence regions bed file", RESOURCE_DIR);
        configBuilder.addPrefixedPath(PANEL_BED, false, "Panel regions bed file", RESOURCE_DIR);
        configBuilder.addPrefixedPath(HOTSPOTS, false, "Hotspots", RESOURCE_DIR);
        DriverGenePanelConfig.addGenePanelOption(configBuilder, false);
        configBuilder.addFlag(PANEL_ONLY, "Only examine panel for variants");

        configBuilder.addFlag(RUN_TINC, "Run TINC routine");
        TincConfig.registerCommonConfig(configBuilder);

        addEnsemblDir(configBuilder);
    }

    public SageCallConfig()
    {
        Common = new SageConfig(false);
        TumorIds = Lists.newArrayList();
        TumorBams = Lists.newArrayList();
        HighConfidenceBed = "highConf";
        PanelBed = "panel";
        Hotspots = "hotspots";
        PanelOnly = false;
        RunTinc = false;
        mResourceDir = "";
    }
}
