package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SampleAnalyserConfig
{
    // reference data
    public final String RefSampleDataFile;
    public final String RefSnvCountsFile;
    public final String RefSvPercFile;
    public final String RefSigContribData;
    public final String RefFeaturePrevFile;
    public final String RefTraitPercFile;
    public final String RefTraitRateFile;
    public final String RefSnvPosFreqFile;
    public final String RefFeatureAvgFile;
    public final String RefRnaExpFile;

    // sample data, if not sourced from the database
    public final String SampleDataFile;
    public final String SampleFeatureFile;
    public final String SampleTraitsFile;
    public final String SampleSnvCountsFile;
    public final String SampleSnvPosFreqFile;
    public final String SampleSigContribFile;
    public final String SampleSvFile;
    public final String SampleRnaExpFile;

    public final List<CategoryType> Categories; // to run, all if empty

    // database access
    public final DatabaseAccess DbAccess;

    public final String OutputDir;
    public final String OutputFileId;

    // config strings
    public static final String CATEGORIES = "categories";
    public static final String SPECIFIC_SAMPLE_DATA = "sample_data";
    public static final String SAMPLE_DATA_FILE = "sample_data_file";
    private static final String SAMPLE_FEAT_FILE = "sample_feature_file";
    private static final String SAMPLE_TRAITS_FILE = "sample_traits_file";
    private static final String SAMPLE_SNV_COUNTS_FILE = "sample_snv_counts_file";
    private static final String SAMPLE_SNV_POS_FREQ_FILE = "sample_snv_pos_freq_file";
    private static final String SAMPLE_SIG_CONTRIB_FILE = "sample_sig_contrib_file";
    private static final String SAMPLE_SV_FILE = "sample_sv_file";
    private static final String SAMPLE_RNA_EXP_FILE = "sample_rna_exp_file";

    public static final String REF_SAMPLE_DATA_FILE = "ref_sample_data_file";
    public static final String REF_SNV_COUNTS_FILE = "ref_snv_counts_file";
    private static final String REF_SIG_CONTRIB_FILE = "ref_sig_contrib_file";
    private static final String REF_FEAT_PREV_FILE = "ref_feature_prev_file";
    private static final String REF_TRAIT_PERC_FILE = "ref_trait_perc_file";
    private static final String REF_TRAIT_RATE_FILE = "ref_trait_rate_file";
    private static final String REF_SV_PERC_FILE = "ref_sv_perc_file";
    private static final String REF_SNV_POS_FREQ_FILE = "ref_snv_pos_freq_file";
    private static final String REF_FEAT_AVG_FILE = "ref_feature_avg_file";
    private static final String REF_RNA_EXP_FILE = "ref_rna_exp_file";

    public static final String OUTPUT_FILE_ID = "output_file_id";
    public static final String LOG_DEBUG = "log_debug";

    public static final Logger CUP_LOGGER = LogManager.getLogger(SampleAnalyserConfig.class);

    public static final String CANCER_SUBTYPE_OTHER = "Other";
    public static final String DATA_DELIM = ",";
    public static final String SUBSET_DELIM = ";";

    public SampleAnalyserConfig(final CommandLine cmd)
    {
        Categories = Lists.newArrayList();

        if(cmd.hasOption(CATEGORIES))
        {
            Categories.addAll(Arrays.stream(cmd.getOptionValue(CATEGORIES)
                    .split(SUBSET_DELIM, -1))
                    .map(x -> CategoryType.valueOf(x))
                    .collect(Collectors.toList()));
        }

        SampleDataFile = cmd.getOptionValue(SAMPLE_DATA_FILE, "");
        SampleTraitsFile = cmd.getOptionValue(SAMPLE_TRAITS_FILE, "");
        SampleFeatureFile = cmd.getOptionValue(SAMPLE_FEAT_FILE, "");
        SampleSnvCountsFile = cmd.getOptionValue(SAMPLE_SNV_COUNTS_FILE, "");
        SampleSnvPosFreqFile = cmd.getOptionValue(SAMPLE_SNV_POS_FREQ_FILE, "");
        SampleSigContribFile = cmd.getOptionValue(SAMPLE_SIG_CONTRIB_FILE, "");
        SampleSvFile = cmd.getOptionValue(SAMPLE_SV_FILE, "");
        SampleRnaExpFile = cmd.getOptionValue(SAMPLE_RNA_EXP_FILE, "");

        RefSampleDataFile = cmd.getOptionValue(REF_SAMPLE_DATA_FILE, "");
        RefSnvCountsFile = cmd.getOptionValue(REF_SNV_COUNTS_FILE, "");
        RefSigContribData = cmd.getOptionValue(REF_SIG_CONTRIB_FILE, "");
        RefFeaturePrevFile = cmd.getOptionValue(REF_FEAT_PREV_FILE, "");
        RefTraitPercFile = cmd.getOptionValue(REF_TRAIT_PERC_FILE, "");
        RefSvPercFile = cmd.getOptionValue(REF_SV_PERC_FILE, "");
        RefTraitRateFile = cmd.getOptionValue(REF_TRAIT_RATE_FILE, "");
        RefSnvPosFreqFile = cmd.getOptionValue(REF_SNV_POS_FREQ_FILE, "");
        RefFeatureAvgFile = cmd.getOptionValue(REF_FEAT_AVG_FILE, "");
        RefRnaExpFile = cmd.getOptionValue(REF_RNA_EXP_FILE, "");

        OutputDir = parseOutputDir(cmd);
        OutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID, "");

        DbAccess = createDatabaseAccess(cmd);
    }

    public boolean isValid()
    {
        return !OutputDir.isEmpty();
    }

    public boolean runCategory(final CategoryType type)
    {
        return Categories.isEmpty() || Categories.contains(type);
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
        options.addOption(CATEGORIES, true, "Optional: list of categories to run classifications, separated by ';' (default=all)");
        options.addOption(SPECIFIC_SAMPLE_DATA, true, "Specific sample in form 'SampleId;CancerType;CancerSubtype' (last 2 optional)");
        options.addOption(SAMPLE_DATA_FILE, true, "Sample data file");
        options.addOption(SAMPLE_SNV_COUNTS_FILE, true, "Sample SNV counts");
        options.addOption(SAMPLE_SNV_POS_FREQ_FILE, true, "Sample SNV position frequence counts");
        options.addOption(SAMPLE_FEAT_FILE, true, "Sample drivers");
        options.addOption(SAMPLE_TRAITS_FILE, true, "Sample traits");
        options.addOption(SAMPLE_SV_FILE, true, "Sample SV data");
        options.addOption(SAMPLE_SIG_CONTRIB_FILE, true, "Sample signature contributions");
        options.addOption(SAMPLE_RNA_EXP_FILE, true, "Sample RNA gene expression TPMs");

        options.addOption(REF_SAMPLE_DATA_FILE, true, "Reference sample data");
        options.addOption(REF_SNV_COUNTS_FILE, true, "Reference SNV sample counts");
        options.addOption(REF_SIG_CONTRIB_FILE, true, "SNV signatures");
        options.addOption(REF_FEAT_PREV_FILE, true, "Reference driver prevalence");
        options.addOption(REF_SV_PERC_FILE, true, "Reference SV percentiles file");
        options.addOption(REF_TRAIT_PERC_FILE, true, "Reference traits percentiles file");
        options.addOption(REF_TRAIT_RATE_FILE, true, "Reference traits rates file");
        options.addOption(REF_SNV_POS_FREQ_FILE, true, "Reference SNV position frequency file");
        options.addOption(REF_FEAT_AVG_FILE, true, "Reference features per sample file");
        options.addOption(REF_RNA_EXP_FILE, true, "Reference RNA gene expression file");

        addDatabaseCmdLineArgs(options);

        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(OUTPUT_FILE_ID, true, "Output file ID");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
    }

}
