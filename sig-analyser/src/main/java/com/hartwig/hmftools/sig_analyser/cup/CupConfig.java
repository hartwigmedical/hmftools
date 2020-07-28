package com.hartwig.hmftools.sig_analyser.cup;

import static com.hartwig.hmftools.sig_analyser.SigAnalyser.LOG_DEBUG;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.OUTPUT_DIR;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.OUTPUT_FILE_ID;

import java.io.File;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CupConfig
{
    public final String SnvSampleCountsFile;
    public final String SampleDataFile;
    public final String SnvSignaturesFile;
    public final String OutputDir;
    public final String OutputFileId;

    private static final String SNV_SAMPLE_COUNTS = "snv_sample_counts";
    private static final String SNV_SIGNATURES = "snv_signatures";
    private static final String SAMPLE_DATA_FILE = "sample_data_file";

    public static final Logger CUP_LOGGER = LogManager.getLogger(CupConfig.class);

    public static final String CANCER_SUBTYPE_OTHER = "Other";

    public CupConfig(final CommandLine cmd)
    {
        SampleDataFile = cmd.getOptionValue(SAMPLE_DATA_FILE, "");
        SnvSampleCountsFile = cmd.getOptionValue(SNV_SAMPLE_COUNTS, "");
        SnvSignaturesFile = cmd.getOptionValue(SNV_SIGNATURES, "");

        String outputDir = cmd.getOptionValue(OUTPUT_DIR);
        if (!outputDir.endsWith(File.separator))
            outputDir += File.separator;

        OutputDir = outputDir;
        OutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID, "");
    }

    public boolean isValid()
    {
        return !SampleDataFile.isEmpty() && !OutputDir.isEmpty();
    }

    public String formOutputFilename(final String fileId)
    {
        if(!OutputFileId.isEmpty())
            return OutputDir + fileId + "." + OutputFileId + ".csv";
        else
            return OutputDir + fileId + ".csv";
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE_DATA_FILE, true, "Sample data file");
        options.addOption(SNV_SAMPLE_COUNTS, true, "SNV sample counts");
        options.addOption(SNV_SIGNATURES, true, "SNV signatures");

        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(OUTPUT_FILE_ID, true, "Output file ID");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
    }

}
