package com.hartwig.hmftools.ctdna;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class PvConfig
{
    public final List<String> SampleIds;
    public final String LinxDir;
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

    public static final Logger PV_LOGGER = LogManager.getLogger(PvConfig.class);

    // config strings
    public static final String SAMPLE = "sample";
    private static final String LINX_DIR = "linx_dir";
    private static final String PURPLE_DIR = "purple_dir";
    private static final String REFERENCE_VARIANTS_FILE = "ref_variants";
    private static final String VAF_THRESHOLD = "vaf_min";
    private static final String FRAG_COUNT_THRESHOLD = "frag_count_min";
    private static final String PROBE_COUNT = "probe_count";
    private static final String PROBE_LENGTH = "probe_length";
    private static final String NON_REPORTABLE_SV_COUNT = "non_reportable_sv_count";
    private static final String SUBCLONAL_COUNT = "subclonal_count";
    private static final String WRITE_ALL = "write_all";

    private static final int DEFAULT_PROBE_COUNT = 500;
    private static final int DEFAULT_PROBE_LENGTH = 120;
    private static final double DEFAULT_VAF_MIN = 0.05;
    private static final int DEFAULT_FRAG_COUNT_MIN = 11;
    private static final int DEFAULT_NON_REPORTABLE_SV_COUNT = 30;
    public static final double DEFAULT_GC_THRESHOLD_MIN = 0.3;
    public static final double DEFAULT_GC_THRESHOLD_MAX = 0.6;
    public static final double DEFAULT_MAPPABILITY_MIN = 0.5;
    public static final double DEFAULT_REPEAT_COUNT_MAX = 3;
    public static final double DEFAULT_SUBCLONAL_LIKELIHOOD_MIN = 0.95;

    public static final int MAX_INSERT_BASES = 60;

    public PvConfig(final CommandLine cmd)
    {
        SampleIds = Lists.newArrayList();

        if(cmd.hasOption(SAMPLE))
        {
            SampleIds.add(cmd.getOptionValue(SAMPLE));
        }
        else
        {
            SampleIds.addAll(loadSampleIdsFile(cmd));
        }

        PurpleDir = cmd.getOptionValue(PURPLE_DIR);
        LinxDir = cmd.getOptionValue(LINX_DIR);
        RefGenomeFile = cmd.getOptionValue(REF_GENOME);
        OutputDir = parseOutputDir(cmd);
        ReferenceVariantsFile = cmd.getOptionValue(REFERENCE_VARIANTS_FILE);

        RefGenVersion = cmd.hasOption(REF_GENOME_VERSION) ? RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION)) : V37;

        PV_LOGGER.info("refGenome({}), purpleDir({}) linxDir({})", RefGenVersion, PurpleDir, LinxDir);
        PV_LOGGER.info("output({})", OutputDir);

        ProbeCount = Integer.parseInt(cmd.getOptionValue(PROBE_COUNT, String.valueOf(DEFAULT_PROBE_COUNT)));
        ProbeLength = Integer.parseInt(cmd.getOptionValue(PROBE_LENGTH, String.valueOf(DEFAULT_PROBE_LENGTH)));
        NonReportableSvCount = Integer.parseInt(cmd.getOptionValue(NON_REPORTABLE_SV_COUNT, String.valueOf(DEFAULT_NON_REPORTABLE_SV_COUNT)));
        SubclonalCount = Integer.parseInt(cmd.getOptionValue(SUBCLONAL_COUNT, "0"));
        VafMin = Double.parseDouble(cmd.getOptionValue(VAF_THRESHOLD, String.valueOf(DEFAULT_VAF_MIN)));
        FragmentCountMin = Integer.parseInt(cmd.getOptionValue(FRAG_COUNT_THRESHOLD, String.valueOf(DEFAULT_FRAG_COUNT_MIN)));
        WriteAll = cmd.hasOption(WRITE_ALL);
        Threads = parseThreads(cmd);
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
            if(!checkFilePath(PURPLE_DIR, PurpleDir, true))
                return false;

            if(!checkFilePath(LINX_DIR, LinxDir, true))
                return false;
        }

        if(!checkFilePath(REF_GENOME, RefGenomeFile, true))
            return false;

        if(!checkFilePath(REFERENCE_VARIANTS_FILE, ReferenceVariantsFile, false))
            return false;

        if(SampleIds.isEmpty())
        {
            PV_LOGGER.error("missing sampleId config");
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
                PV_LOGGER.error("missing config: {}", config);
                return false;
            }

            return true;
        }

        if(!Files.exists(Paths.get(filePath)))
        {
            PV_LOGGER.error("invalid {} path: {}", config, filePath);
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
        RefGenomeFile = "";
        OutputDir = "";
        ReferenceVariantsFile = "";
        RefGenVersion = V37;
        WriteAll = false;
        Threads = 1;
    }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();

        options.addOption(SAMPLE, true, "Tumor sample ID");
        addSampleIdFile(options);
        options.addOption(LINX_DIR, true, "Linx directory");
        options.addOption(PURPLE_DIR, true, "Purple directory");
        addRefGenomeConfig(options);
        options.addOption(REFERENCE_VARIANTS_FILE, true, "Reference variants file");
        options.addOption(VAF_THRESHOLD, true, "VAF threshold, default: " + DEFAULT_VAF_MIN);
        options.addOption(FRAG_COUNT_THRESHOLD, true, "Fragment count threshold, default: " + DEFAULT_FRAG_COUNT_MIN);
        options.addOption(PROBE_COUNT, true, "Probe count, default: " + DEFAULT_PROBE_COUNT);
        options.addOption(PROBE_LENGTH, true, "Probe length, default: " + DEFAULT_PROBE_LENGTH);
        options.addOption(NON_REPORTABLE_SV_COUNT, true, "Max count of non-reportable SVs, default: " + DEFAULT_NON_REPORTABLE_SV_COUNT);
        options.addOption(SUBCLONAL_COUNT, true, "Max count of subclonal mutations, default(0)");
        options.addOption(WRITE_ALL, false, "Write all variants to file");

        addOutputOptions(options);
        addLoggingOptions(options);
        addThreadOptions(options, 1);

        return options;
    }
}
