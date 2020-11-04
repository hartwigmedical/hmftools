package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_DRIVER_AVG;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_FEATURE_PREV;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_GENE_EXP_CANCER;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_GENE_EXP_PERC;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_GENE_EXP_SAMPLE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_POS_FREQ_COUNTS;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SAMPLE_DATA;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SIG_PERC;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SNV_COUNTS;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SV_PERC;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_TRAIT_PERC;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_TRAIT_RATES;
import static com.hartwig.hmftools.cup.common.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.common.CategoryType.GENE_EXP;
import static com.hartwig.hmftools.cup.common.CategoryType.SAMPLE_TRAIT;
import static com.hartwig.hmftools.cup.common.CategoryType.SNV;
import static com.hartwig.hmftools.cup.common.CategoryType.SV;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.rna.RnaExpression;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CuppaConfig
{
    public final List<CategoryType> IncludedCategories;

    // reference data
    public final String RefDataDir;

    public final String RefSampleDataFile;
    public final String RefSnvCountsFile;
    public final String RefSvPercFile;
    public final String RefSigContributionFile;
    public final String RefFeaturePrevFile;
    public final String RefTraitPercFile;
    public final String RefTraitRateFile;
    public final String RefSnvPosFreqFile;
    public final String RefDriverAvgFile;

    public final String RefGeneExpCancerFile;
    public final String RefGeneExpSampleFile;
    public final String RefGeneExpPercFile;

    // sample data, if not sourced from the database
    public final String SampleDataDir;

    public final String SampleDataFile;
    public final String SampleFeatureFile;
    public final String SampleTraitsFile;
    public final String SampleSnvCountsFile;
    public final String SampleSnvPosFreqFile;
    public final String SampleSigContribFile;
    public final String SampleSomaticVcf;
    public final String SampleSvFile;
    public final String SampleRnaExpFile;

    // database access
    public final DatabaseAccess DbAccess;

    public final boolean WriteSimilarities;
    public final String OutputDir;
    public final String OutputFileId;

    // config strings
    public static final String CATEGORIES = "categories";

    public static final String SAMPLE_DATA_DIR = "sample_data_dir";

    public static final String SPECIFIC_SAMPLE_DATA = "sample_data";
    public static final String SAMPLE_DATA_FILE = "sample_data_file";
    private static final String SAMPLE_FEAT_FILE = "sample_feature_file";
    private static final String SAMPLE_TRAITS_FILE = "sample_traits_file";
    private static final String SAMPLE_SNV_COUNTS_FILE = "sample_snv_counts_file";
    private static final String SAMPLE_SNV_POS_FREQ_FILE = "sample_snv_pos_freq_file";
    private static final String SAMPLE_SIG_CONTRIB_FILE = "sample_sig_contrib_file";
    private static final String SAMPLE_SOMATIC_VCF = "sample_somatic_vcf";
    private static final String SAMPLE_SV_FILE = "sample_sv_file";
    private static final String SAMPLE_RNA_EXP_FILE = "sample_gene_exp_file";

    public static final String REF_DATA_DIR = "ref_data_dir";

    public static final String REF_SAMPLE_DATA_FILE = "ref_sample_data_file";
    public static final String REF_SNV_COUNTS_FILE = "ref_snv_counts_file";
    private static final String REF_SNV_POS_FREQ_FILE = "ref_snv_pos_freq_file";
    private static final String REF_SIG_CONTRIB_FILE = "ref_sig_contrib_file";
    private static final String REF_FEAT_PREV_FILE = "ref_feature_prev_file";
    private static final String REF_DRIVER_AVG_FILE = "ref_feature_avg_file";
    private static final String REF_TRAIT_PERC_FILE = "ref_trait_perc_file";
    private static final String REF_TRAIT_RATE_FILE = "ref_trait_rate_file";
    private static final String REF_SV_PERC_FILE = "ref_sv_perc_file";
    private static final String REF_RNA_GENE_EXP_CANCER_FILE = "ref_gene_exp_cancer_file";
    private static final String REF_RNA_GENE_EXP_SAMPLE_FILE = "ref_gene_exp_sample_file";
    private static final String REF_RNA_GENE_EXP_CANCER_PERC_FILE = "ref_gene_exp_cancer_perc_file";

    public static final String WRITE_SIMS = "write_similarities";
    public static final String OUTPUT_FILE_ID = "output_file_id";
    public static final String LOG_DEBUG = "log_debug";

    public static final Logger CUP_LOGGER = LogManager.getLogger(CuppaConfig.class);

    public static final String CANCER_SUBTYPE_OTHER = "Other";
    public static final String DATA_DELIM = ",";
    public static final String SUBSET_DELIM = ";";

    public CuppaConfig(final CommandLine cmd)
    {
        IncludedCategories = Lists.newArrayList();

        if(cmd.hasOption(CATEGORIES))
        {
            if(cmd.getOptionValue(CATEGORIES).equals("ALL"))
            {
                // will run all classifiers
            }
            else
            {
                final String[] categories = cmd.getOptionValue(CATEGORIES).split(";");
                Arrays.stream(categories).forEach(x -> IncludedCategories.add(CategoryType.valueOf(x)));
            }
        }
        else
        {
            IncludedCategories.add(SNV);
            IncludedCategories.add(SV);
            IncludedCategories.add(FEATURE);
            IncludedCategories.add(SAMPLE_TRAIT);
        }

        CUP_LOGGER.info("running classifiers: {}", IncludedCategories.toString());

        SampleDataDir = checkAddDirSeparator(cmd.getOptionValue(SAMPLE_DATA_DIR, ""));

        SampleDataFile = cmd.getOptionValue(SAMPLE_DATA_FILE, "");
        SampleTraitsFile = cmd.getOptionValue(SAMPLE_TRAITS_FILE, "");
        SampleFeatureFile = cmd.getOptionValue(SAMPLE_FEAT_FILE, "");
        SampleSnvCountsFile = cmd.getOptionValue(SAMPLE_SNV_COUNTS_FILE, "");
        SampleSnvPosFreqFile = cmd.getOptionValue(SAMPLE_SNV_POS_FREQ_FILE, "");
        SampleSigContribFile = cmd.getOptionValue(SAMPLE_SIG_CONTRIB_FILE, "");
        SampleSvFile = cmd.getOptionValue(SAMPLE_SV_FILE, "");
        SampleRnaExpFile = cmd.getOptionValue(SAMPLE_RNA_EXP_FILE, "");
        SampleSomaticVcf = cmd.getOptionValue(SAMPLE_SOMATIC_VCF, "");

        RefDataDir = checkAddDirSeparator(cmd.getOptionValue(REF_DATA_DIR, ""));

        RefSampleDataFile = getRefDataFile(cmd, REF_SAMPLE_DATA_FILE, REF_FILE_SAMPLE_DATA);
        RefSnvCountsFile = getRefDataFile(cmd, REF_SNV_COUNTS_FILE, REF_FILE_SNV_COUNTS);
        RefSigContributionFile = getRefDataFile(cmd, REF_SIG_CONTRIB_FILE, REF_FILE_SIG_PERC);
        RefFeaturePrevFile = getRefDataFile(cmd, REF_FEAT_PREV_FILE, REF_FILE_FEATURE_PREV);
        RefTraitPercFile = getRefDataFile(cmd, REF_TRAIT_PERC_FILE, REF_FILE_TRAIT_PERC);
        RefTraitRateFile = getRefDataFile(cmd, REF_TRAIT_RATE_FILE, REF_FILE_TRAIT_RATES);
        RefSvPercFile = getRefDataFile(cmd, REF_SV_PERC_FILE, REF_FILE_SV_PERC);
        RefSnvPosFreqFile = getRefDataFile(cmd, REF_SNV_POS_FREQ_FILE, REF_FILE_POS_FREQ_COUNTS);
        RefDriverAvgFile = getRefDataFile(cmd, REF_DRIVER_AVG_FILE, REF_FILE_DRIVER_AVG);

        RefGeneExpCancerFile = getRefDataFile(cmd, REF_RNA_GENE_EXP_CANCER_FILE, REF_FILE_GENE_EXP_CANCER);
        RefGeneExpSampleFile = getRefDataFile(cmd, REF_RNA_GENE_EXP_SAMPLE_FILE, REF_FILE_GENE_EXP_SAMPLE);
        RefGeneExpPercFile = getRefDataFile(cmd, REF_RNA_GENE_EXP_CANCER_PERC_FILE, REF_FILE_GENE_EXP_PERC);

        OutputDir = parseOutputDir(cmd);
        OutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID, "");

        WriteSimilarities = Boolean.parseBoolean(cmd.getOptionValue(WRITE_SIMS, "true"));

        DbAccess = createDatabaseAccess(cmd);
    }

    private String getRefDataFile(final CommandLine cmd, final String configStr, final String defaultFilename)
    {
        final String fileName = cmd.getOptionValue(configStr, defaultFilename);
        return RefDataDir + fileName;
    }

    public boolean isValid()
    {
        return !OutputDir.isEmpty();
    }

    public boolean runClassifier(final CategoryType type) { return IncludedCategories.isEmpty() || IncludedCategories.contains(type); }

    public String formOutputFilename(final String fileId)
    {
        String outputFile = OutputDir + "CUP";

        if(!OutputFileId.isEmpty())
            outputFile += "." + OutputFileId;

        return outputFile + "." + fileId + ".csv";
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(
                CATEGORIES, true,
                "Categories to run analysis on, separated by ';' from: SNV, SV, SAMPLE_TRAIT, GENE_EXP and FEATURE");

        options.addOption(SPECIFIC_SAMPLE_DATA, true, "Specific sample in form 'SampleId;CancerType;CancerSubtype' (last 2 optional)");
        options.addOption(SAMPLE_DATA_DIR, true, "Directory containing standard sample files from pipeline");

        options.addOption(SAMPLE_DATA_FILE, true, "Sample data file");

        options.addOption(SAMPLE_SNV_COUNTS_FILE, true, "Sample SNV counts");
        options.addOption(SAMPLE_SNV_POS_FREQ_FILE, true, "Sample SNV position frequence counts");
        options.addOption(SAMPLE_SIG_CONTRIB_FILE, true, "Sample signature contributions");

        options.addOption(SAMPLE_FEAT_FILE, true, "Cohort features file (drivers, fusions and viruses)");
        options.addOption(SAMPLE_TRAITS_FILE, true, "Cohort sample traits file");
        options.addOption(SAMPLE_SV_FILE, true, "Cohort SV data");
        options.addOption(SAMPLE_RNA_EXP_FILE, true, "Sample RNA gene expression TPMs");
        options.addOption(SAMPLE_SOMATIC_VCF, true, "Sample somatic VCF");

        options.addOption(REF_DATA_DIR, true, "Reference data directory");

        options.addOption(REF_SAMPLE_DATA_FILE, true, "Reference sample data, default: " + REF_FILE_SAMPLE_DATA);
        options.addOption(REF_SNV_COUNTS_FILE, true, "Reference SNV sample counts, default: " + REF_FILE_SNV_COUNTS);
        options.addOption(REF_SIG_CONTRIB_FILE, true, "SNV signatures, default: " + REF_FILE_SIG_PERC);
        options.addOption(REF_FEAT_PREV_FILE, true, "Reference driver prevalence, default: " + REF_FILE_FEATURE_PREV);
        options.addOption(REF_SV_PERC_FILE, true, "Reference SV percentiles file, default: " + REF_FILE_SV_PERC);
        options.addOption(REF_TRAIT_PERC_FILE, true, "Reference traits percentiles file, default: " + REF_FILE_TRAIT_PERC);
        options.addOption(REF_TRAIT_RATE_FILE, true, "Reference traits rates file, default: " + REF_FILE_TRAIT_RATES);
        options.addOption(REF_SNV_POS_FREQ_FILE, true, "Reference SNV position frequency file, default: " + REF_FILE_POS_FREQ_COUNTS);
        options.addOption(REF_DRIVER_AVG_FILE, true, "Reference features per sample file, default: " + REF_FILE_DRIVER_AVG);
        options.addOption(REF_RNA_GENE_EXP_CANCER_FILE, true, "Reference RNA cancer gene expression file, default: " + REF_FILE_GENE_EXP_CANCER);
        options.addOption(REF_RNA_GENE_EXP_SAMPLE_FILE, true, "Reference RNA sample gene expression file, default: " + REF_FILE_GENE_EXP_SAMPLE);
        options.addOption(REF_RNA_GENE_EXP_CANCER_PERC_FILE, true, "Reference RNA cancer percentiles gene expression file, default: " + REF_FILE_GENE_EXP_PERC);

        options.addOption(WRITE_SIMS, true, "Cohort-only - write top-20 CSS similarities to file");

        addDatabaseCmdLineArgs(options);
        RnaExpression.addCmdLineArgs(options);

        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(OUTPUT_FILE_ID, true, "Output file ID");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
    }

}
