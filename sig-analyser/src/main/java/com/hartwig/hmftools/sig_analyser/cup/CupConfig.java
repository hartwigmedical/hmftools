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
    public final String DriverPrevFile;
    public final String SampleDriversFile;
    public final String OutputDir;
    public final String OutputFileId;

    private static final String SNV_SAMPLE_COUNTS = "snv_sample_counts";
    private static final String SNV_SIGNATURES = "snv_signatures";
    private static final String SAMPLE_DATA_FILE = "sample_data_file";
    private static final String DRIVER_PREV_FILE = "driver_prev_file";
    private static final String SAMPLE_DRIVERS_FILE = "sample_drivers_file";

    public static final Logger CUP_LOGGER = LogManager.getLogger(CupConfig.class);

    public static final String CANCER_SUBTYPE_OTHER = "Other";
    public static final String DATA_DELIM = ",";

    public CupConfig(final CommandLine cmd)
    {
        SampleDataFile = cmd.getOptionValue(SAMPLE_DATA_FILE, "");
        SnvSampleCountsFile = cmd.getOptionValue(SNV_SAMPLE_COUNTS, "");
        SnvSignaturesFile = cmd.getOptionValue(SNV_SIGNATURES, "");
        DriverPrevFile = cmd.getOptionValue(DRIVER_PREV_FILE, "");
        SampleDriversFile = cmd.getOptionValue(SAMPLE_DRIVERS_FILE, "");

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
        String outputFile = OutputDir + "CUP";

        if(!OutputFileId.isEmpty())
            outputFile += "." + OutputFileId;

        return outputFile + "." + fileId + ".csv";
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE_DATA_FILE, true, "Sample data file");
        options.addOption(SNV_SAMPLE_COUNTS, true, "SNV sample counts");
        options.addOption(SNV_SIGNATURES, true, "SNV signatures");
        options.addOption(DRIVER_PREV_FILE, true, "Driver prevalence");
        options.addOption(SAMPLE_DRIVERS_FILE, true, "Sample drivers");

        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(OUTPUT_FILE_ID, true, "Output file ID");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
    }

}
