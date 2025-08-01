package com.hartwig.hmftools.panelbuilder.wisp;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.panelbuilder.wisp.ProbeConstants.DEFAULT_FRAG_COUNT_MIN;
import static com.hartwig.hmftools.panelbuilder.wisp.ProbeConstants.DEFAULT_FRAG_COUNT_MIN_LOWER;
import static com.hartwig.hmftools.panelbuilder.wisp.ProbeConstants.DEFAULT_GC_THRESHOLD_MAX;
import static com.hartwig.hmftools.panelbuilder.wisp.ProbeConstants.DEFAULT_GC_THRESHOLD_MAX_LOWER;
import static com.hartwig.hmftools.panelbuilder.wisp.ProbeConstants.DEFAULT_GC_THRESHOLD_MIN;
import static com.hartwig.hmftools.panelbuilder.wisp.ProbeConstants.DEFAULT_GC_THRESHOLD_MIN_LOWER;
import static com.hartwig.hmftools.panelbuilder.wisp.ProbeConstants.DEFAULT_PROBE_COUNT;
import static com.hartwig.hmftools.panelbuilder.wisp.ProbeConstants.DEFAULT_VAF_MIN;

import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ProbeConfig
{
    public final String SampleId;
    public final String LinxDir;
    public final String LinxGermlineDir;
    public final String PurpleDir;

    public final String ReferenceVariantsFile;

    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final double VafMin;
    public final int FragmentCountMin;
    public final double GcRatioLimitMin;
    public final double GcRatioLimitMax;

    // for second-round considerations
    public final int FragmentCountMinLower;
    public final double GcRatioLimitLowerMin;
    public final double GcRatioLimitLowerMax;

    public final int ProbeCount;
    public final int NonReportableSvCount;
    public final int SubclonalCount;
    public final boolean WriteAll;

    public final String OutputDir;
    public final String OutputId;

    // config strings
    private static final String REFERENCE_VARIANTS_FILE = "ref_variants";
    private static final String VAF_THRESHOLD = "vaf_min";
    private static final String FRAG_COUNT_THRESHOLD = "frag_count_min";
    private static final String FRAG_COUNT_THRESHOLD_LOWER = "frag_count_min_lower";
    private static final String PROBE_COUNT = "probe_count";
    private static final String NON_REPORTABLE_SV_COUNT = "non_reportable_sv_count";
    private static final String SUBCLONAL_COUNT = "subclonal_count";
    private static final String WRITE_ALL = "write_all";

    private static final Logger LOGGER = LogManager.getLogger(ProbeConfig.class);

    public ProbeConfig(final ConfigBuilder configBuilder)
    {
        SampleId = configBuilder.getValue(SAMPLE);
        PurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        LinxDir = configBuilder.getValue(LINX_DIR_CFG);
        LinxGermlineDir = configBuilder.getValue(LINX_GERMLINE_DIR_CFG);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        ReferenceVariantsFile = configBuilder.getValue(REFERENCE_VARIANTS_FILE);

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        LOGGER.info("refGenome({}), purpleDir({}) linxDir({})", RefGenVersion, PurpleDir, LinxDir);
        LOGGER.info("output({})", OutputDir);

        ProbeCount = configBuilder.getInteger(PROBE_COUNT);
        NonReportableSvCount = configBuilder.getInteger(NON_REPORTABLE_SV_COUNT);
        SubclonalCount = configBuilder.getInteger(SUBCLONAL_COUNT);
        VafMin = configBuilder.getDecimal(VAF_THRESHOLD);
        FragmentCountMin = configBuilder.getInteger(FRAG_COUNT_THRESHOLD);
        FragmentCountMinLower = configBuilder.getInteger(FRAG_COUNT_THRESHOLD_LOWER);

        GcRatioLimitMin = DEFAULT_GC_THRESHOLD_MIN;
        GcRatioLimitMax = DEFAULT_GC_THRESHOLD_MAX;
        GcRatioLimitLowerMin = DEFAULT_GC_THRESHOLD_MIN_LOWER;
        GcRatioLimitLowerMax = DEFAULT_GC_THRESHOLD_MAX_LOWER;

        WriteAll = configBuilder.hasFlag(WRITE_ALL);
    }

    public static String getSampleFilePath(final String sampleId, final String filePath)
    {
        return convertWildcardSamplePath(filePath, sampleId);
    }

    public void checkSampleDirectories(final String sampleId)
    {
        String purpleDir = getSampleFilePath(sampleId, PurpleDir);

        // allow Linx inputs to be optional
        String linxDir = getSampleFilePath(sampleId, LinxDir);
        String linxGermlineDir = getSampleFilePath(sampleId, LinxGermlineDir);

        if(!Files.exists(Paths.get(purpleDir))
                || (linxDir != null && !Files.exists(Paths.get(linxDir)))
                || (linxGermlineDir != null && !Files.exists(Paths.get(linxGermlineDir))))
        {
            String error = format("sample(%s) missing Purple or Linx directories", sampleId);
            LOGGER.error(error);
            throw new RuntimeException(error);
        }
    }

    public ProbeConfig(int probeCount, int nonReportableSvCount, double vafMin, int fragmentCountMin)
    {
        ProbeCount = probeCount;
        NonReportableSvCount = nonReportableSvCount;
        SubclonalCount = 0;
        VafMin = vafMin;
        FragmentCountMin = fragmentCountMin;
        FragmentCountMinLower = DEFAULT_FRAG_COUNT_MIN_LOWER;
        GcRatioLimitMin = DEFAULT_GC_THRESHOLD_MIN;
        GcRatioLimitMax = DEFAULT_GC_THRESHOLD_MAX;
        GcRatioLimitLowerMin = DEFAULT_GC_THRESHOLD_MIN_LOWER;
        GcRatioLimitLowerMax = DEFAULT_GC_THRESHOLD_MAX_LOWER;
        SampleId = "";
        PurpleDir = "";
        LinxDir = "";
        LinxGermlineDir = "";
        RefGenomeFile = "";
        OutputDir = "";
        OutputId = "";
        ReferenceVariantsFile = "";
        RefGenVersion = V37;
        WriteAll = false;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        configBuilder.addPath(LINX_DIR_CFG, false, LINX_DIR_DESC);
        configBuilder.addPath(LINX_GERMLINE_DIR_CFG, false, LINX_GERMLINE_DIR_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        addRefGenomeConfig(configBuilder, true);
        configBuilder.addPath(REFERENCE_VARIANTS_FILE, false, "Reference variants file");

        configBuilder.addDecimal(VAF_THRESHOLD, "VAF threshold", DEFAULT_VAF_MIN);

        configBuilder.addInteger(FRAG_COUNT_THRESHOLD, "Fragment count threshold", DEFAULT_FRAG_COUNT_MIN);

        configBuilder.addInteger(
                FRAG_COUNT_THRESHOLD_LOWER, "Fragment count lower threshold for second-round selection",
                DEFAULT_FRAG_COUNT_MIN_LOWER);

        configBuilder.addInteger(PROBE_COUNT, "Probe count", DEFAULT_PROBE_COUNT);

        configBuilder.addInteger(NON_REPORTABLE_SV_COUNT, "Max count of non-reportable SVs", 0);
        configBuilder.addInteger(SUBCLONAL_COUNT, "Max count of subclonal mutations", 0);

        configBuilder.addFlag(WRITE_ALL, "Write all variants to file");

        addOutputOptions(configBuilder);
    }
}
