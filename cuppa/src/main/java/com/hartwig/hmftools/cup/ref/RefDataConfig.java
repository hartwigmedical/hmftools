package com.hartwig.hmftools.cup.ref;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.REF_SAMPLE_DATA_FILE;
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
    public final String RefSampleDataFile;
    public final String RefSampleTraitsFile;
    public final String RefSigContribsFile;

    // config strings
    public static final String REF_SAMPLE_TRAITS_FILE = "ref_sample_traits_file";
    public static final String REF_SIG_CONTRIBS_FILE = "ref_sig_contribs_file";

    public static final String DB_FILE_DELIM = "\t";

    public RefDataConfig(final CommandLine cmd)
    {
        RefSampleDataFile = cmd.getOptionValue(REF_SAMPLE_DATA_FILE, "");
        RefSampleTraitsFile = cmd.getOptionValue(REF_SAMPLE_TRAITS_FILE, "");
        RefSigContribsFile = cmd.getOptionValue(REF_SIG_CONTRIBS_FILE, "");

        String outputDir = cmd.getOptionValue(OUTPUT_DIR);
        if (!outputDir.endsWith(File.separator))
            outputDir += File.separator;

        OutputDir = outputDir;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(REF_SAMPLE_DATA_FILE, true, "Sample ref data file");
        options.addOption(REF_SAMPLE_TRAITS_FILE, true, "Ref sample traits data file");
        options.addOption(REF_SIG_CONTRIBS_FILE, true, "Ref signature contirbutions data file");

        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
    }

}
