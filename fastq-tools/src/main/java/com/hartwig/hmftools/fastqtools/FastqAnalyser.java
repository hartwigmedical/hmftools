package com.hartwig.hmftools.fastqtools;

import static java.lang.String.format;

import static com.hartwig.hmftools.fastqtools.FastqCommon.FQ_LOGGER;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;

public class FastqAnalyser
{
    private final String mFastqFile;
    private final String mOutputDir;

    // config
    private static final String FASTQ_FILE = "fastq_file";

    public FastqAnalyser(final CommandLine cmd)
    {
        mFastqFile = cmd.getOptionValue(FASTQ_FILE);
        mOutputDir = parseOutputDir(cmd);
    }

    public void run()
    {
        if(mOutputDir == null || mFastqFile == null)
            System.exit(1);

        FQ_LOGGER.info("Starting Fastq Analyser with file: {}", mFastqFile);


        long startTimeMs = System.currentTimeMillis();


        long timeTakenMs = System.currentTimeMillis() - startTimeMs;
        double timeTakeMins = timeTakenMs / 60000.0;

        FQ_LOGGER.info("Fastq analysis complete, mins({})", format("%.3f", timeTakeMins));
    }

    public static void main(@NotNull final String[] args)
    {
        // final VersionInfo version = new VersionInfo("fastq-tools.version");
        // FQ_LOGGER.info("BamTools version: {}", version.version());

        final Options options = new Options();

        addOutputOptions(options);
        addLoggingOptions(options);
        //addThreadOptions(options);
        //addRefGenomeConfig(options);;
        options.addOption(FASTQ_FILE, true, "Fastq file path");

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            FastqAnalyser fastqAnalyser = new FastqAnalyser(cmd);
            fastqAnalyser.run();
        }
        catch(ParseException e)
        {
            FQ_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("FastqAnalyser", options);
            System.exit(1);
        }
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
