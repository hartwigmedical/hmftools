package com.hartwig.hmftools.ctdna.probe;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.DEFAULT_PROBE_LENGTH;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class PvConfig
{
    public final List<String> SampleIds;
    public final String LinxDir;
    public final String LinxGermlineDir;
    public final String PurpleDir;

    public final String ReferenceVariantsFile;

    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final double VafMin;
    public final int FragmentCountMin;
    public final int ProbeCount;
    public final int ProbeLength;
    public final int NonReportableSvCount;
    public final int SubclonalCount;
    public final boolean WriteAll;

    public final int Threads;
    public final String OutputDir;
    public final String OutputId;

    // config strings
    public static final String SAMPLE = "sample";
    private static final String LINX_GERMLINE_DIR = "linx_germline_dir";
    private static final String REFERENCE_VARIANTS_FILE = "ref_variants";
    private static final String VAF_THRESHOLD = "vaf_min";
    private static final String FRAG_COUNT_THRESHOLD = "frag_count_min";
    private static final String PROBE_COUNT = "probe_count";
    private static final String PROBE_LENGTH = "probe_length";
    private static final String NON_REPORTABLE_SV_COUNT = "non_reportable_sv_count";
    private static final String SUBCLONAL_COUNT = "subclonal_count";
    private static final String WRITE_ALL = "write_all";

    private static final int DEFAULT_PROBE_COUNT = 500;
    private static final double DEFAULT_VAF_MIN = 0.05;
    private static final int DEFAULT_FRAG_COUNT_MIN = 11;
    private static final int DEFAULT_NON_REPORTABLE_SV_COUNT = 30;
    public static final double DEFAULT_GC_THRESHOLD_MIN = 0.3;
    public static final double DEFAULT_GC_THRESHOLD_MAX = 0.6;
    public static final double DEFAULT_MAPPABILITY_MIN = 0.5;
    public static final double DEFAULT_REPEAT_COUNT_MAX = 3;
    public static final double DEFAULT_SUBCLONAL_LIKELIHOOD_MIN = 0.95;
    public static final int DEFAULT_SV_BREAKENDS_PER_GENE = 5;

    public static final int MAX_INSERT_BASES = 60;
    public static final int MAX_INDEL_LENGTH = 32;

    public PvConfig(final ConfigBuilder configBuilder)
    {
        SampleIds = Lists.newArrayList();

        if(configBuilder.hasValue(SAMPLE))
        {
            SampleIds.add(configBuilder.getValue(SAMPLE));
        }
        else
        {
            SampleIds.addAll(loadSampleIdsFile(configBuilder));
        }

        PurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        LinxDir = configBuilder.getValue(LINX_DIR_CFG);
        LinxGermlineDir = configBuilder.getValue(LINX_GERMLINE_DIR);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        ReferenceVariantsFile = configBuilder.getValue(REFERENCE_VARIANTS_FILE);

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        CT_LOGGER.info("refGenome({}), purpleDir({}) linxDir({})", RefGenVersion, PurpleDir, LinxDir);
        CT_LOGGER.info("output({})", OutputDir);

        ProbeCount = configBuilder.getInteger(PROBE_COUNT);
        ProbeLength = configBuilder.getInteger(PROBE_LENGTH);
        NonReportableSvCount = configBuilder.getInteger(NON_REPORTABLE_SV_COUNT);
        SubclonalCount = configBuilder.getInteger(SUBCLONAL_COUNT);
        VafMin = configBuilder.getDecimal(VAF_THRESHOLD);
        FragmentCountMin = configBuilder.getInteger(FRAG_COUNT_THRESHOLD);
        WriteAll = configBuilder.hasFlag(WRITE_ALL);
        Threads = parseThreads(configBuilder);
    }

    public boolean isMultiSample() { return SampleIds.size() > 1; }
    public String sample() { return SampleIds.get(0); }

    public static String getSampleFilePath(final String sampleId, final String filePath)
    {
        return filePath.replaceAll("\\*", sampleId);
    }

    public boolean isValid()
    {
        if(SampleIds.size() == 1)
        {
            if(!checkFilePath(PURPLE_DIR_CFG, PurpleDir, true))
                return false;

            if(!checkFilePath(LINX_DIR_CFG, LinxDir, true) || !checkFilePath(LINX_GERMLINE_DIR, LinxDir, true))
                return false;
        }

        if(!checkFilePath(REF_GENOME, RefGenomeFile, true))
            return false;

        if(!checkFilePath(REFERENCE_VARIANTS_FILE, ReferenceVariantsFile, false))
            return false;

        if(SampleIds.isEmpty())
        {
            CT_LOGGER.error("missing sampleId config");
            return false;
        }

        return true;
    }

    private static boolean checkFilePath(final String config, final String filePath, boolean required)
    {
        if(filePath == null)
        {
            if(required)
            {
                CT_LOGGER.error("missing config: {}", config);
                return false;
            }

            return true;
        }

        if(!Files.exists(Paths.get(filePath)))
        {
            CT_LOGGER.error("invalid {} path: {}", config, filePath);
            return false;
        }

        return true;
    }

    public PvConfig(int probeLength, int probeCount, int nonReportableSvCount, double vafMin, int fragmentCountMin)
    {
        ProbeCount = probeCount;
        ProbeLength = probeLength;
        NonReportableSvCount = nonReportableSvCount;
        SubclonalCount = 0;
        VafMin = vafMin;
        FragmentCountMin = fragmentCountMin;
        SampleIds = Lists.newArrayList();
        PurpleDir = "";
        LinxDir = "";
        LinxGermlineDir = "";
        RefGenomeFile = "";
        OutputDir = "";
        OutputId = "";
        ReferenceVariantsFile = "";
        RefGenVersion = V37;
        WriteAll = false;
        Threads = 1;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, false, "Tumor sample ID");
        addSampleIdFile(configBuilder, false);
        configBuilder.addConfigItem(LINX_DIR_CFG, true, LINX_DIR_DESC);
        configBuilder.addConfigItem(LINX_GERMLINE_DIR, true, "Linx germline directory");
        configBuilder.addConfigItem(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        addRefGenomeConfig(configBuilder, true);
        configBuilder.addConfigItem(REFERENCE_VARIANTS_FILE, false, "Reference variants file");

        configBuilder.addDecimal(VAF_THRESHOLD, "VAF threshold", DEFAULT_VAF_MIN);

        configBuilder.addInteger(
                FRAG_COUNT_THRESHOLD, "Fragment count threshold",
                DEFAULT_FRAG_COUNT_MIN);

        configBuilder.addInteger(PROBE_COUNT, "Probe count", DEFAULT_PROBE_COUNT);
        configBuilder.addInteger(PROBE_LENGTH, "Probe length", DEFAULT_PROBE_LENGTH);

        configBuilder.addInteger(
                NON_REPORTABLE_SV_COUNT,"Max count of non-reportable SVs", DEFAULT_NON_REPORTABLE_SV_COUNT);

        configBuilder.addInteger(SUBCLONAL_COUNT, "Max count of subclonal mutations, default(0)", 0);

        configBuilder.addFlag(WRITE_ALL, "Write all variants to file");

        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
    }
}
