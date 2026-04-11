package com.hartwig.hmftools.linx.cohort;

import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigItem.enumValueSelectionAsStr;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class CohortConfig
{
    public final String InputDir;
    public final String OutputDir;
    public final String OutputId;

    public final String PonFile;
    public final int PonMargin;

    public final List<WriteType> WriteTypes;

    private static final String INPUT_DIR = "input_dir";
    private static final String WRITE_TYPES = "write_types";
    private static final String PON_FILE = "pon_file";
    private static final String PON_MARGIN = "pon_margin";

    private static final int DEFAULT_PON_MARGIN = 10;

    public CohortConfig(final ConfigBuilder configBuilder)
    {
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        InputDir = checkAddDirSeparator(configBuilder.getValue(INPUT_DIR, OutputDir));

        PonFile = configBuilder.getValue(PON_FILE);
        PonMargin = configBuilder.getInteger(PON_MARGIN);

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
        configBuilder.addPath(PON_FILE, false, "Pre-built PON file");

        configBuilder.addConfigItem(WRITE_TYPES, false, enumValueSelectionAsStr(WriteType.values(), "Write types"));
        configBuilder.addInteger(PON_MARGIN, "PON margin for matching", DEFAULT_PON_MARGIN);

        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
