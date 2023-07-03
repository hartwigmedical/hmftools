package com.hartwig.hmftools.cobalt.norm;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.addGcProfilePath;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.AMBER_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.AMBER_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.COBALT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.COBALT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_COLUMN;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class NormalisationConfig
{
    public final List<String> SampleIds;
    public final String AmberDir;
    public final String CobaltWgsDir;
    public final String CobaltPanelDir; // unnormalised files
    public final String TargetRegionsBed;
    public final String GcProfile;
    public final String OutputFile;
    public final String DetailedFile;
    public final RefGenomeVersion RefGenVersion;

    private final Map<String,String> mPanelToWgsSampleIdMappings; // if required, mapping from panel to WGS

    private static final String COBALT_WGS_DIR = "cobalt_wgs_dir";
    private static final String TARGET_REGIONS_BED = "target_regions_bed";
    private static final String OUTPUT_FILE = "output_file";
    private static final String DETAILED_OUTPUT = "detailed_file";

    private static final String WGS_SAMPLE_ID = "WgsSampleId";

    public NormalisationConfig(final ConfigBuilder configBuilder)
    {
        SampleIds = Lists.newArrayList();
        mPanelToWgsSampleIdMappings = Maps.newHashMap();
        loadSampleIds(configBuilder);

        CobaltPanelDir = configBuilder.getValue(COBALT_DIR_CFG);
        CobaltWgsDir = configBuilder.getValue(COBALT_WGS_DIR);
        AmberDir = configBuilder.getValue(AMBER_DIR_CFG);
        GcProfile = configBuilder.getValue(GC_PROFILE);
        TargetRegionsBed = configBuilder.getValue(TARGET_REGIONS_BED);
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        OutputFile = configBuilder.getValue(OUTPUT_FILE);
        DetailedFile = configBuilder.getValue(DETAILED_OUTPUT);
    }

    public String getWgsSampleId(final String sampleId)
    {
        return mPanelToWgsSampleIdMappings.containsKey(sampleId) ? mPanelToWgsSampleIdMappings.get(sampleId) : sampleId;
    }

    private void loadSampleIds(final ConfigBuilder configBuilder)
    {
        if(!configBuilder.hasValue(SAMPLE_ID_FILE))
            return;

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(configBuilder.getValue(SAMPLE_ID_FILE)));
            String header = lines.get(0);
            lines.remove(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, CSV_DELIM);

            int sampleIndex = fieldsIndexMap.get(SAMPLE_ID_COLUMN);
            Integer wgsSampleIndex = fieldsIndexMap.get(WGS_SAMPLE_ID);

            for(String line : lines)
            {
                String[] values = line.split(CSV_DELIM, -1);

                String sampleId = values[sampleIndex];

                if(sampleId.startsWith("#"))
                    continue;

                SampleIds.add(sampleId);

                if(wgsSampleIndex != null)
                {
                    mPanelToWgsSampleIdMappings.put(sampleId, values[wgsSampleIndex]);
                }
            }

            CB_LOGGER.info("loaded {} samples from file", SampleIds.size());
        }
        catch(IOException e)
        {
            CB_LOGGER.error("failed to load sample IDs: {}", e.toString());
        }
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        addSampleIdFile(configBuilder, true);
        configBuilder.addPath(AMBER_DIR_CFG, true, AMBER_DIR_DESC);
        configBuilder.addPath(COBALT_DIR_CFG, true, COBALT_DIR_DESC);
        configBuilder.addPath(COBALT_WGS_DIR, false, "Path to cobalt WGS files");
        configBuilder.addConfigItem(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        configBuilder.addPath(TARGET_REGIONS_BED, true, "Target regions BED file");
        configBuilder.addRequiredConfigItem(OUTPUT_FILE, "Output normalisation file");
        configBuilder.addFlag(DETAILED_OUTPUT, "Output normalisation file");
        addGcProfilePath(configBuilder, true);
        addLoggingOptions(configBuilder);
    }
}
