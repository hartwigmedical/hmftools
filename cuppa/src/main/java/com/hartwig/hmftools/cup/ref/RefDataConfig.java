package com.hartwig.hmftools.cup.ref;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.cup.CuppaConfig.LOG_DEBUG;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SNV_COUNTS_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SNV_SAMPLE_POS_FREQ_FILE;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.cup.rna.RefGeneExpression;
import com.hartwig.hmftools.cup.somatics.RefSomatics;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

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
    public final String RefFeaturesFile;
    public final String RefGeneExpFile;
    public final String RefAltSjFile;
    public final String RefFeatureOverrideFile;

    public final DatabaseAccess DbAccess;

    public final boolean WriteCohortFiles; // re-write data sourced from database or flat files into single cohort files

    // config strings

    // cohort files to avoid database data extraction
    public static final String REF_SAMPLE_TRAITS_FILE = "ref_sample_traits_file";
    public static final String REF_SIG_CONTRIBS_FILE = "ref_sig_contribs_file";
    public static final String REF_SV_DATA_FILE = "ref_sv_data_file";
    public static final String REF_FEATURES_FILE = "ref_features_file";

    public static final String REF_GENE_EXP_DATA_FILE = "ref_gene_exp_file";
    public static final String REF_ALT_SJ_DATA_FILE = "ref_alt_sj_file";
    public static final String REF_FEATURE_OVERRIDE_FILE = "feature_override_file";
    public static final String GENDER_RATES = "gender_rates";
    public static final String WRITE_COHORT_FILES = "write_cohort_files";

    public static final String FILE_DELIM = ",";

    public RefDataConfig(final CommandLine cmd)
    {
        RefSampleDataFile = cmd.getOptionValue(REF_SAMPLE_DATA_FILE, "");
        RefSampleTraitsFile = cmd.getOptionValue(REF_SAMPLE_TRAITS_FILE, "");
        RefSigContribsFile = cmd.getOptionValue(REF_SIG_CONTRIBS_FILE, "");
        RefSampleSvDataFile = cmd.getOptionValue(REF_SV_DATA_FILE, "");
        RefSnvPositionDataFile = cmd.getOptionValue(REF_SNV_SAMPLE_POS_FREQ_FILE, "");
        RefSnvCountsFile = cmd.getOptionValue(REF_SNV_COUNTS_FILE, "");
        RefFeaturesFile = cmd.getOptionValue(REF_FEATURES_FILE, "");
        RefGeneExpFile = cmd.getOptionValue(REF_GENE_EXP_DATA_FILE, "");
        RefAltSjFile = cmd.getOptionValue(REF_ALT_SJ_DATA_FILE, "");
        RefFeatureOverrideFile = cmd.getOptionValue(REF_FEATURE_OVERRIDE_FILE, "");

        DbAccess = createDatabaseAccess(cmd);

        OutputDir = parseOutputDir(cmd);

        WriteCohortFiles = cmd.hasOption(WRITE_COHORT_FILES);
    }

    public static final List<String> parseFileSet(final String filenames)
    {
        return Arrays.stream(filenames.split(FILE_DELIM, -1)).collect(Collectors.toList());
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(REF_SAMPLE_DATA_FILE, true, "Ref sample data file");
        options.addOption(REF_SAMPLE_TRAITS_FILE, true, "Ref sample traits data file");
        options.addOption(REF_SIG_CONTRIBS_FILE, true, "Ref signature contributions data file");
        options.addOption(REF_SV_DATA_FILE, true, "Ref sample SV data file");
        options.addOption(REF_SNV_SAMPLE_POS_FREQ_FILE, true, "Ref SNV position frequency matrix data file");
        options.addOption(REF_SNV_COUNTS_FILE, true, "Ref SNV trinucleotide matrix data file");
        options.addOption(REF_FEATURES_FILE, true, "Ref sample features data file");
        options.addOption(REF_GENE_EXP_DATA_FILE, true, "Ref sample RNA gene expression cohort data file");
        options.addOption(REF_ALT_SJ_DATA_FILE, true, "Ref sample RNA alternate splice junction cohort data file");
        options.addOption(REF_FEATURE_OVERRIDE_FILE, true, "Ref feature override data file");
        options.addOption(GENDER_RATES, true, "Gender-rate overrides - format CancerType;MalePerc;FemalePerc, etc");
        options.addOption(WRITE_COHORT_FILES, false, "Re-write ref data as cohort files");

        addDatabaseCmdLineArgs(options);

        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");

        RefGeneExpression.addCmdLineArgs(options);
        RefSomatics.addCmdLineArgs(options);
    }

}
