package com.hartwig.hmftools.cup.ref;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.LOG_DEBUG;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.REF_SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.REF_SNV_COUNTS_FILE;

import java.io.File;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class RefDataConfig
{
    public final String OutputDir;

    // reference data
    public final String RefSampleDataFile;
    public final String RefSampleTraitsFile;
    public final String RefSampleSvDataFile;
    public final String RefSigContribsFile;
    public final String RefSnvPositionDataFile;
    public final String RefSnvCountsFile;

    // config strings
    public static final String REF_SAMPLE_TRAITS_FILE = "ref_sample_traits_file";
    public static final String REF_SIG_CONTRIBS_FILE = "ref_sig_contribs_file";
    public static final String REF_SV_DATA_FILE = "ref_sv_data_file";
    public static final String REF_SNV_POS_DATA_FILE = "ref_snv_pos_file";

    public static final String DB_FILE_DELIM = "\t";

    public RefDataConfig(final CommandLine cmd)
    {
        RefSampleDataFile = cmd.getOptionValue(REF_SAMPLE_DATA_FILE, "");
        RefSampleTraitsFile = cmd.getOptionValue(REF_SAMPLE_TRAITS_FILE, "");
        RefSigContribsFile = cmd.getOptionValue(REF_SIG_CONTRIBS_FILE, "");
        RefSampleSvDataFile = cmd.getOptionValue(REF_SV_DATA_FILE, "");
        RefSnvPositionDataFile = cmd.getOptionValue(REF_SNV_POS_DATA_FILE, "");
        RefSnvCountsFile = cmd.getOptionValue(REF_SNV_COUNTS_FILE, "");

        OutputDir = parseOutputDir(cmd);
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(REF_SAMPLE_DATA_FILE, true, "Ref sample data file");
        options.addOption(REF_SAMPLE_TRAITS_FILE, true, "Ref sample traits data file");
        options.addOption(REF_SIG_CONTRIBS_FILE, true, "Ref signature contirbutions data file");
        options.addOption(REF_SV_DATA_FILE, true, "Ref sample SV data file");
        options.addOption(REF_SNV_POS_DATA_FILE, true, "Ref sample SNV position bucket data file");
        options.addOption(REF_SNV_COUNTS_FILE, true, "Ref sample SNV position bucket data file");

        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
    }

}
