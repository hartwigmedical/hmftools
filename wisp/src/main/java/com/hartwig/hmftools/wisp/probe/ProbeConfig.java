package com.hartwig.hmftools.wisp.probe;

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
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.common.CommonUtils.DEFAULT_PROBE_LENGTH;
import static com.hartwig.hmftools.wisp.probe.ProbeConstants.DEFAULT_FRAG_COUNT_MIN;
import static com.hartwig.hmftools.wisp.probe.ProbeConstants.DEFAULT_FRAG_COUNT_MIN_LOWER;
import static com.hartwig.hmftools.wisp.probe.ProbeConstants.DEFAULT_GC_THRESHOLD_MAX;
import static com.hartwig.hmftools.wisp.probe.ProbeConstants.DEFAULT_GC_THRESHOLD_MAX_LOWER;
import static com.hartwig.hmftools.wisp.probe.ProbeConstants.DEFAULT_GC_THRESHOLD_MIN;
import static com.hartwig.hmftools.wisp.probe.ProbeConstants.DEFAULT_GC_THRESHOLD_MIN_LOWER;
import static com.hartwig.hmftools.wisp.probe.ProbeConstants.DEFAULT_PROBE_COUNT;
import static com.hartwig.hmftools.wisp.probe.ProbeConstants.DEFAULT_VAF_MIN;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class ProbeConfig
{
    public final List<String> SampleIds;
    public final Map<String,List<String>> BatchSampleIds;
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
    public final int ProbeLength;
    public final int NonReportableSvCount;
    public final int SubclonalCount;
    public final boolean WriteAll;
    public final boolean AllowMissing;

    public final int Threads;
    public final String OutputDir;
    public final String OutputId;

    // config strings
    private static final String REFERENCE_VARIANTS_FILE = "ref_variants";
    private static final String VAF_THRESHOLD = "vaf_min";
    private static final String FRAG_COUNT_THRESHOLD = "frag_count_min";
    private static final String FRAG_COUNT_THRESHOLD_LOWER = "frag_count_min_lower";
    private static final String PROBE_COUNT = "probe_count";
    private static final String PROBE_LENGTH = "probe_length";
    private static final String NON_REPORTABLE_SV_COUNT = "non_reportable_sv_count";
    private static final String SUBCLONAL_COUNT = "subclonal_count";
    private static final String WRITE_ALL = "write_all";
    private static final String ALLOW_MISSING = "allow_missing";
    private static final String SAMPLE_BATCH_COUNT = "sample_batch_count";

    public static final String NO_BATCH_ID = "NONE";

    public ProbeConfig(final ConfigBuilder configBuilder)
    {
        SampleIds = Lists.newArrayList();
        BatchSampleIds = Maps.newHashMap();

        if(configBuilder.hasValue(SAMPLE))
        {
            SampleIds.add(configBuilder.getValue(SAMPLE));
        }
        else
        {
            loadSampleIdsFile(configBuilder.getValue(SAMPLE_ID_FILE));
        }

        PurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        LinxDir = configBuilder.getValue(LINX_DIR_CFG);
        LinxGermlineDir = configBuilder.getValue(LINX_GERMLINE_DIR_CFG);
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
        FragmentCountMinLower = configBuilder.getInteger(FRAG_COUNT_THRESHOLD_LOWER);

        GcRatioLimitMin = DEFAULT_GC_THRESHOLD_MIN;
        GcRatioLimitMax = DEFAULT_GC_THRESHOLD_MAX;
        GcRatioLimitLowerMin = DEFAULT_GC_THRESHOLD_MIN_LOWER;
        GcRatioLimitLowerMax = DEFAULT_GC_THRESHOLD_MAX_LOWER;

        WriteAll = configBuilder.hasFlag(WRITE_ALL);
        AllowMissing = configBuilder.hasFlag(ALLOW_MISSING);
        Threads = parseThreads(configBuilder);
    }

    public boolean isMultiSample() { return SampleIds.size() > 1; }
    public String sample() { return SampleIds.get(0); }

    public static String getSampleFilePath(final String sampleId, final String filePath)
    {
        return convertWildcardSamplePath(filePath, sampleId);
    }

    public boolean isValid()
    {
        if(SampleIds.isEmpty())
        {
            CT_LOGGER.error("missing sampleId config");
            return false;
        }

        return true;
    }

    public boolean checkSampleDirectories(final String sampleId)
    {
        String purpleDir = getSampleFilePath(sampleId, PurpleDir);

        // allow Linx inputs to be optional
        String linxDir = getSampleFilePath(sampleId, LinxDir);
        String linxGermlineDir = getSampleFilePath(sampleId, LinxGermlineDir);

        if(!Files.exists(Paths.get(purpleDir))
        || (linxDir != null && !Files.exists(Paths.get(linxDir)))
        || (linxGermlineDir != null && !Files.exists(Paths.get(linxGermlineDir))))
        {
            if(!AllowMissing)
            {
                CT_LOGGER.error("sample({}) missing Purple or Linx directories", sampleId);
                System.exit(1);
            }

            CT_LOGGER.warn("sample({}) missing Purple or Linx directories", sampleId);
            return false;
        }

        return true;
    }

    private void loadSampleIdsFile(final String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(new File(filename).toPath());
            String header = lines.get(0);
            lines.remove(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, CSV_DELIM);

            int idIndex = fieldsIndexMap.get(FLD_SAMPLE_ID);
            Integer batchIndex = fieldsIndexMap.get("BatchId");

            for(String line : lines)
            {
                String[] values = line.split(CSV_DELIM, -1);

                String sampleId = values[idIndex];

                String batchId = batchIndex != null ? values[batchIndex] : NO_BATCH_ID;

                List<String> sampleIds = BatchSampleIds.get(batchId);

                if(sampleIds == null)
                {
                    sampleIds = Lists.newArrayList();
                    BatchSampleIds.put(batchId, sampleIds);
                }

                sampleIds.add(sampleId);

                SampleIds.add(sampleId);
            }

            if(!BatchSampleIds.isEmpty())
            {
                CT_LOGGER.info("loaded {} samples in {} batches from file({})", SampleIds.size(), BatchSampleIds.size(), filename);
            }
            else
            {
                CT_LOGGER.info("loaded {} samples from file({})", SampleIds.size(), filename);
            }
        }
        catch (IOException e)
        {
            CT_LOGGER.error("failed to read sample ID file({}): {}", filename, e.toString());
        }
    }

    public ProbeConfig(int probeLength, int probeCount, int nonReportableSvCount, double vafMin, int fragmentCountMin)
    {
        ProbeCount = probeCount;
        ProbeLength = probeLength;
        NonReportableSvCount = nonReportableSvCount;
        SubclonalCount = 0;
        VafMin = vafMin;
        FragmentCountMin = fragmentCountMin;
        FragmentCountMinLower = DEFAULT_FRAG_COUNT_MIN_LOWER;
        GcRatioLimitMin = DEFAULT_GC_THRESHOLD_MIN;
        GcRatioLimitMax = DEFAULT_GC_THRESHOLD_MAX;
        GcRatioLimitLowerMin = DEFAULT_GC_THRESHOLD_MIN_LOWER;
        GcRatioLimitLowerMax = DEFAULT_GC_THRESHOLD_MAX_LOWER;
        SampleIds = Lists.newArrayList();
        BatchSampleIds = Maps.newHashMap();
        PurpleDir = "";
        LinxDir = "";
        LinxGermlineDir = "";
        RefGenomeFile = "";
        OutputDir = "";
        OutputId = "";
        ReferenceVariantsFile = "";
        RefGenVersion = V37;
        WriteAll = false;
        AllowMissing = false;
        Threads = 1;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        addSampleIdFile(configBuilder, false);
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
        configBuilder.addInteger(PROBE_LENGTH, "Probe length", DEFAULT_PROBE_LENGTH);
        configBuilder.addInteger(SAMPLE_BATCH_COUNT, "Sample batching count", 1);

        configBuilder.addInteger(NON_REPORTABLE_SV_COUNT,"Max count of non-reportable SVs", 0);
        configBuilder.addInteger(SUBCLONAL_COUNT, "Max count of subclonal mutations", 0);

        configBuilder.addFlag(WRITE_ALL, "Write all variants to file");
        configBuilder.addFlag(ALLOW_MISSING, "Continue on missing sample data");

        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
    }
}
