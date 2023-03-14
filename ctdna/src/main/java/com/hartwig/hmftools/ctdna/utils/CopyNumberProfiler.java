package com.hartwig.hmftools.ctdna.utils;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;

import com.hartwig.hmftools.ctdna.purity.CopyNumberProfile;
import com.hartwig.hmftools.ctdna.purity.PurityConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class CopyNumberProfiler
{
    private static final String SAMPLE = "sample";

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Sample ID");
        PurityConfig.addCommandLineOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        PurityConfig purityConfig = new PurityConfig(cmd);
        String sampleId = cmd.getOptionValue(SAMPLE);

        CopyNumberProfile cnProfile = new CopyNumberProfile(purityConfig);
        cnProfile.loadSampleData(sampleId);

        CT_LOGGER.info("Sample VCF analyser complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
