package com.hartwig.hmftools.esvee.depth;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificRegions;
import static com.hartwig.hmftools.common.bam.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.region.UnmappedRegions.UNMAP_REGIONS_FILE;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.INPUT_VCF;
import static com.hartwig.hmftools.esvee.common.FileCommon.INPUT_VCF_DESC;
import static com.hartwig.hmftools.esvee.common.FileCommon.RAW_VCF_SUFFIX;
import static com.hartwig.hmftools.esvee.common.FileCommon.formEsveeInputFilename;
import static com.hartwig.hmftools.esvee.common.FileCommon.parseBamFiles;
import static com.hartwig.hmftools.esvee.common.FileCommon.parseSampleIds;
import static com.hartwig.hmftools.esvee.prep.PrepConfig.BAM_FILE;
import static com.hartwig.hmftools.esvee.prep.PrepConfig.BAM_FILE_DESC;
import static com.hartwig.hmftools.esvee.prep.PrepConfig.SAMPLE_ID_DESC;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.BamUtils;
import com.hartwig.hmftools.common.region.UnmappedRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.apache.commons.cli.ParseException;

import htsjdk.samtools.ValidationStringency;

public class DepthConfig
{
    public final String InputVcf;
    public final String OutputDir;
    public final String OutputId;
    public final List<String> SampleIds;
    public final List<String> BamFiles;
    public final String RefGenome;
    public final String UnmapRegionsFile;
    public final RefGenomeVersion RefGenVersion;
    public final double VafCap;
    public final int ProximityDistance;
    public final String VcfTagPrefix;
    public final ValidationStringency BamStringency;

    public final int Threads;
    public final double PerfLogTime;
    public final List<ChrBaseRegion> SpecificRegions;

    private static final String PROXIMITY_DISTANCE = "proximity_distance";
    private static final String VAF_CAP = "vaf_cap";
    private static final String PERF_LOG_TIME = "perf_log_time";
    public static final String VCF_TAG_PREFIX = "vcf_tag_prefix";

    protected static final int DEFAULT_PROXIMITY_DISTANCE = 2000;
    protected static final double DEFAULT_VAF_CAP = 0.001;

    public DepthConfig(final ConfigBuilder configBuilder)
    {
        SampleIds = parseSampleIds(configBuilder);
        BamFiles = parseBamFiles(configBuilder);

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        if(configBuilder.hasValue(INPUT_VCF))
            InputVcf = configBuilder.getValue(INPUT_VCF);
        else
            InputVcf = formEsveeInputFilename(OutputDir, SampleIds.get(0), RAW_VCF_SUFFIX, OutputId);

        RefGenome = configBuilder.getValue(REF_GENOME);
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        VafCap = configBuilder.getDecimal(VAF_CAP);
        ProximityDistance = configBuilder.getInteger(PROXIMITY_DISTANCE);
        BamStringency = BamUtils.validationStringency(configBuilder);
        PerfLogTime = configBuilder.getDecimal(PERF_LOG_TIME);

        UnmapRegionsFile = configBuilder.getValue(UNMAP_REGIONS_FILE);

        Threads = parseThreads(configBuilder);

        SpecificRegions = Lists.newArrayList();

        try
        {
            SpecificRegions.addAll(loadSpecificRegions(configBuilder));
        }
        catch(ParseException e)
        {
            SV_LOGGER.error("failed to load specific regions");
        }

        VcfTagPrefix = configBuilder.getValue(VCF_TAG_PREFIX);
    }

    public String sampleId() { return SampleIds.get(0); }

    public String getVcfTag(final String vcfTag)
    {
        return VcfTagPrefix != null ? format("%s_%s", VcfTagPrefix, vcfTag) : vcfTag;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(INPUT_VCF, false, INPUT_VCF_DESC);
        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_ID_DESC); // not required, see comment in prep config
        configBuilder.addPaths(BAM_FILE, false, BAM_FILE_DESC);
        configBuilder.addConfigItem(VCF_TAG_PREFIX, "VCF tag prefix for testing & comparison");
        addRefGenomeConfig(configBuilder, true);
        UnmappedRegions.registerConfig(configBuilder);

        configBuilder.addDecimal(VAF_CAP, "Ref support depth limit as function of variant fragments", DEFAULT_VAF_CAP);
        configBuilder.addInteger(PROXIMITY_DISTANCE, "Proximity distance to group variants", DEFAULT_PROXIMITY_DISTANCE);

        addValidationStringencyOption(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);
        addThreadOptions(configBuilder);
        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.addDecimal(PERF_LOG_TIME, "Performance log time threshold (seconds)", 0);
    }

    public DepthConfig(double vcfCap, int proximityDistance)
    {
        InputVcf = "";
        OutputDir = "";
        OutputId = "";

        SampleIds = Lists.newArrayList();
        BamFiles = Lists.newArrayList();

        RefGenome = "";
        RefGenVersion = V37;
        VafCap = vcfCap;
        ProximityDistance = proximityDistance;
        UnmapRegionsFile = null;
        BamStringency = ValidationStringency.STRICT;
        PerfLogTime = 0;
        Threads = 0;

        SpecificRegions = Lists.newArrayList();
        VcfTagPrefix = "";
    }
}
