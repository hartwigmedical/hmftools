package com.hartwig.hmftools.linx.cohort;

import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigItem.enumValueSelectionAsStr;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;

import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class CohortConfig
{
    public final String InputDir;
    public final String OutputDir;
    public final String OutputId;

    public final String PonFile;
    public final int PonMargin;
    public final boolean RestrictNonPon;
    public final boolean WritePonType;

    public final List<String> SampleIds;

    public final List<WriteType> WriteTypes;

    private static final String INPUT_DIR = "input_dir";
    private static final String WRITE_TYPES = "write_types";
    private static final String PON_FILE = "pon_file";
    private static final String PON_MARGIN = "pon_margin";
    private static final String RESTRICT_NON_PON  = "restrict_non_pon";
    private static final String WRITE_PON_TYPE  = "write_pon_type";

    private static final int DEFAULT_PON_MARGIN = 10;

    public CohortConfig(final ConfigBuilder configBuilder)
    {
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        InputDir = checkAddDirSeparator(configBuilder.getValue(INPUT_DIR, OutputDir));

        PonFile = configBuilder.getValue(PON_FILE);
        PonMargin = configBuilder.getInteger(PON_MARGIN);
        RestrictNonPon = configBuilder.hasFlag(RESTRICT_NON_PON);
        WritePonType = configBuilder.hasFlag(WRITE_PON_TYPE);

        SampleIds = loadSampleIdsFile(configBuilder);

        if(!SampleIds.isEmpty())
        {
            LNX_LOGGER.info("restricting to {} samples", SampleIds.size());
        }

        WriteTypes = WriteType.fromConfigStr(configBuilder.getValue(WRITE_TYPES));
    }

    public String formOutputFilename(final String fileId)
    {
        String filename = OutputDir + "linx_cohort." + fileId;

        if(OutputId != null)
            filename += "." + OutputId;

        filename += TSV_EXTENSION;
        return filename;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(INPUT_DIR, false, "Directory with Linx cohort files. If blank, uses output directory");
        configBuilder.addConfigItem(PON_FILE, false, "Pre-built PON file. If multiple then name1=file1;name2=file2");

        configBuilder.addConfigItem(WRITE_TYPES, false, enumValueSelectionAsStr(WriteType.values(), "Write types"));
        configBuilder.addInteger(PON_MARGIN, "PON margin for matching", DEFAULT_PON_MARGIN);
        configBuilder.addFlag(RESTRICT_NON_PON, "All output excludes SVs matching the PON");
        configBuilder.addFlag(WRITE_PON_TYPE, "Write PON match type");
        configBuilder.addConfigItem(SAMPLE_ID_FILE, false, SAMPLE_ID_FILE_DESC);

        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
