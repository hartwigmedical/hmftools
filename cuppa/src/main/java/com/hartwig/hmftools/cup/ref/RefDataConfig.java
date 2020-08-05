package com.hartwig.hmftools.cup.ref;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.LOG_DEBUG;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.OUTPUT_DIR;

import java.io.File;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class RefDataConfig
{
    public final String OutputDir;

    // reference data
    public final String SampleDataFile;
    public final String SampleTraitsFile;

    // config strings
    public static final String SAMPLE_TRAITS_FILE = "sample_traits_file";

    public static final String DB_FILE_DELIM = "\t";

    public RefDataConfig(final CommandLine cmd)
    {
        SampleDataFile = cmd.getOptionValue(SAMPLE_DATA_FILE, "");
        SampleTraitsFile = cmd.getOptionValue(SAMPLE_TRAITS_FILE, "");

        String outputDir = cmd.getOptionValue(OUTPUT_DIR);
        if (!outputDir.endsWith(File.separator))
            outputDir += File.separator;

        OutputDir = outputDir;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE_DATA_FILE, true, "Sample ref data file");
        options.addOption(SAMPLE_TRAITS_FILE, true, "Sample traits ref data file");

        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
    }

}
