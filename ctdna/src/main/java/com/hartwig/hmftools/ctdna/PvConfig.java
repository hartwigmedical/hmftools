package com.hartwig.hmftools.ctdna;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;

import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class PvConfig
{
    public final String SampleId;
    public final String LinxDir;
    public final String PurpleDir;

    public final String ActionableVariantsFile;
    public final String PolymorphismsFile;

    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final double VafThreshold;
    public final int FragmentCountThreshold;
    public final int ProbeCount;
    public final int ProbeLength;
    public final int NonReportableSvCount;
    public final boolean WriteAll;

    public final String OutputDir;
    private boolean mIsValid;

    public static final Logger PV_LOGGER = LogManager.getLogger(PvConfig.class);

    // config strings
    public static final String SAMPLE = "sample";
    private static final String LINX_DIR = "linx_dir";
    private static final String PURPLE_DIR = "purple_dir";
    private static final String ACTIONABLE_VARIANTS_FILE = "actionable_variants_file";
    private static final String POLYMORPHISMS_FILE = "polymorphisms_file";
    private static final String VAF_THRESHOLD = "vaf_threshold";
    private static final String FRAG_COUNT_THRESHOLD = "frag_count_threshold";
    private static final String PROBE_COUNT = "probe_count";
    private static final String PROBE_LENGTH = "probe_length";
    private static final String NON_REPORTABLE_SV_COUNT = "non_reportable_sv_count";
    private static final String WRITE_ALL = "write_all";

    private static final int DEFAULT_PROBE_COUNT = 500;
    private static final int DEFAULT_PROBE_LENGTH = 120;
    private static final double DEFAULT_VAF_THRESHOLD = 0.5;
    private static final int DEFAULT_FRAG_COUNT_THRESHOLD = 10;
    private static final int DEFAULT_NON_REPORTABLE_SV_COUNT = 30;

    public static final int MAX_INSERT_BASES = 60;

    public PvConfig(final CommandLine cmd)
    {
        SampleId = cmd.getOptionValue(SAMPLE);
        PurpleDir = cmd.getOptionValue(PURPLE_DIR);
        LinxDir = cmd.getOptionValue(LINX_DIR);
        RefGenomeFile = cmd.getOptionValue(REF_GENOME);
        OutputDir = parseOutputDir(cmd);
        ActionableVariantsFile = cmd.getOptionValue(ACTIONABLE_VARIANTS_FILE);
        PolymorphismsFile = cmd.getOptionValue(POLYMORPHISMS_FILE);

        RefGenVersion = cmd.hasOption(REF_GENOME_VERSION) ? RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION)) : V37;

        PV_LOGGER.info("refGenome({}), purpleDir({}) linxDir({})", RefGenVersion, PurpleDir, LinxDir);
        PV_LOGGER.info("output({})", OutputDir);

        ProbeCount = Integer.parseInt(cmd.getOptionValue(PROBE_COUNT, String.valueOf(DEFAULT_PROBE_COUNT)));
        ProbeLength = Integer.parseInt(cmd.getOptionValue(PROBE_LENGTH, String.valueOf(DEFAULT_PROBE_LENGTH)));
        NonReportableSvCount = Integer.parseInt(cmd.getOptionValue(NON_REPORTABLE_SV_COUNT, String.valueOf(DEFAULT_NON_REPORTABLE_SV_COUNT)));
        VafThreshold = Double.parseDouble(cmd.getOptionValue(VAF_THRESHOLD, String.valueOf(DEFAULT_VAF_THRESHOLD)));
        FragmentCountThreshold = Integer.parseInt(cmd.getOptionValue(FRAG_COUNT_THRESHOLD, String.valueOf(DEFAULT_FRAG_COUNT_THRESHOLD)));
        WriteAll = cmd.hasOption(WRITE_ALL);
    }

    public boolean isValid()
    {
        if(!checkFilePath(PURPLE_DIR, PurpleDir, true))
            return false;

        if(!checkFilePath(LINX_DIR, LinxDir, true))
            return false;

        if(!checkFilePath(REF_GENOME, RefGenomeFile, true))
            return false;

        if(!checkFilePath(ACTIONABLE_VARIANTS_FILE, ActionableVariantsFile, false))
            return false;

        if(!checkFilePath(POLYMORPHISMS_FILE, PolymorphismsFile, false))
            return false;

        if(SampleId == null)
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

    public PvConfig(int probeLength, int probeCount, int nonReportableSvCount, double vafThreshold, int fragmentCountThreshold)
    {
        ProbeCount = probeCount;
        ProbeLength = probeLength;
        NonReportableSvCount = nonReportableSvCount;
        VafThreshold = vafThreshold;
        FragmentCountThreshold = fragmentCountThreshold;
        SampleId = "";
        PurpleDir = "";
        LinxDir = "";
        RefGenomeFile = "";
        OutputDir = "";
        ActionableVariantsFile = "";
        PolymorphismsFile = "";
        RefGenVersion = V37;
        WriteAll = false;
    }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        addOutputOptions(options);
        addLoggingOptions(options);

        options.addOption(SAMPLE, true, "Tumor sample ID");
        options.addOption(LINX_DIR, true, "Linx directory");
        options.addOption(PURPLE_DIR, true, "Purple directory");
        addRefGenomeConfig(options);
        options.addOption(ACTIONABLE_VARIANTS_FILE, true, "Hotspots file");
        options.addOption(POLYMORPHISMS_FILE, true, "Polymorphisms file");
        options.addOption(VAF_THRESHOLD, true, "VAF threshold, default: " + DEFAULT_VAF_THRESHOLD);
        options.addOption(FRAG_COUNT_THRESHOLD, true, "Fragment count threshold, default: " + DEFAULT_FRAG_COUNT_THRESHOLD);
        options.addOption(PROBE_COUNT, true, "Probe count, default: " + DEFAULT_PROBE_COUNT);
        options.addOption(PROBE_LENGTH, true, "Probe length, default: " + DEFAULT_PROBE_LENGTH);
        options.addOption(NON_REPORTABLE_SV_COUNT, true, "Max count of non-reportable SVs, default: " + DEFAULT_NON_REPORTABLE_SV_COUNT);
        options.addOption(WRITE_ALL, false, "Write all variants to file");

        return options;
    }
}
