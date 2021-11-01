package com.hartwig.hmftools.telo.breakend;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.utils.ConfigUtils;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BreakEndAnalyserConfig
{
    public final String SampleId;

    public final String TelbamFile;
    public final String OutputDir;

    public final String SampleType;

    private static final String SAMPLE_ID = "sample_id";
    private static final String TELBAM_FILE = "telbam_file";
    private static final String SAMPLE_TYPE = "sample_type";
    private static final String BOUNDARY_ZONE = "boundary_zone";
    private static final String FUZZY_MATCH_DISTANCE = "fuzzy_match_distance";

    public static final Logger LOGGER = LogManager.getLogger(BreakEndAnalyserConfig.class);

    public BreakEndAnalyserConfig(final CommandLine cmd)
    {
        SampleId = cmd.getOptionValue(SAMPLE_ID);
        TelbamFile = cmd.getOptionValue(TELBAM_FILE);
        OutputDir = parseOutputDir(cmd);
        SampleType = cmd.getOptionValue(SAMPLE_TYPE);
    }

    public boolean isValid()
    {
        if(!checkCreateOutputDir(OutputDir))
        {
            LOGGER.error("failed to create output directory({})", OutputDir);
            return false;
        }

        return true;
    }

    public static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(Option.builder(SAMPLE_ID).hasArg().required().desc("ID of tumor sample").build());
        options.addOption(Option.builder(TELBAM_FILE).hasArg().required().desc("Path to telbam.bam file").build());
        options.addOption(Option.builder(OUTPUT_DIR).hasArg().required().desc("Output directory").build());
        options.addOption(Option.builder(SAMPLE_TYPE).hasArg().required().desc("Type of sample (germline / somatic)").build());
        ConfigUtils.addLoggingOptions(options);

        return options;
    }
}
