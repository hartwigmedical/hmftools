package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.sig_analyser.SigAnalyser.LOG_DEBUG;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.OUTPUT_DIR;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.OUTPUT_FILE_ID;

import java.io.File;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.cup.common.DataType;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SampleAnalyserConfig
{
    public final List<DataType> DataTypes;

    // sample data, if not sourced from the database
    public final String SampleDataFile;
    public final String SampleDriversFile;

    public final String OutputDir;
    public final String OutputFileId;

    // reference data
    public final String SnvSampleCountsFile;
    public final String SnvSignaturesFile;
    public final String DriverPrevFile;

    // config string
    public static final String SPECIFIC_SAMPLE_DATA = "sample_data";
    public static final String SAMPLE_DATA_FILE = "sample_data_file";
    private static final String SAMPLE_DRIVERS_FILE = "sample_drivers_file";

    private static final String SNV_SAMPLE_COUNTS = "snv_sample_counts";
    private static final String SNV_SIGNATURES = "snv_signatures";
    private static final String DRIVER_PREV_FILE = "driver_prev_file";

    public static final Logger CUP_LOGGER = LogManager.getLogger(SampleAnalyserConfig.class);

    public static final String CANCER_SUBTYPE_OTHER = "Other";
    public static final String DATA_DELIM = ",";
    public static final String SUBSET_DELIM = ";";

    public SampleAnalyserConfig(final CommandLine cmd)
    {
        DataTypes = Lists.newArrayList();

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
        return !OutputDir.isEmpty();
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
        options.addOption(SPECIFIC_SAMPLE_DATA, true, "Specific sample in form 'SampleId;CancerType;CancerSubtype' (last 2 optional)");
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
