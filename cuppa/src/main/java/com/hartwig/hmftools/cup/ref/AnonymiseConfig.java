package com.hartwig.hmftools.cup.ref;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.cup.CuppaConfig.LOG_DEBUG;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_DATA_DIR;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_RNA_ALT_SJ_SAMPLE_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_RNA_GENE_EXP_SAMPLE_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SNV_COUNTS_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SNV_SAMPLE_POS_FREQ_FILE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_ALT_SJ_SAMPLE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_GENE_EXP_SAMPLE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SAMPLE_DATA;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SAMPLE_POS_FREQ_COUNTS;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SNV_COUNTS;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class AnonymiseConfig
{
    public final String OutputDir;

    // reference data
    public final String RefDataDir;
    public final String RefSampleDataFile;
    public final String RefSnvSamplePosFreqFile;
    public final String RefSnvCountsFile;
    public final String RefGeneExpSampleFile;
    public final String RefAltSjSampleFile;

    public AnonymiseConfig(final CommandLine cmd)
    {
        RefDataDir = checkAddDirSeparator(cmd.getOptionValue(REF_DATA_DIR, ""));

        RefSampleDataFile = getRefDataFile(cmd, REF_SAMPLE_DATA_FILE, REF_FILE_SAMPLE_DATA);
        RefSnvCountsFile = getRefDataFile(cmd, REF_SNV_COUNTS_FILE, REF_FILE_SNV_COUNTS);
        RefSnvSamplePosFreqFile = getRefDataFile(cmd, REF_SNV_SAMPLE_POS_FREQ_FILE, REF_FILE_SAMPLE_POS_FREQ_COUNTS);
        RefGeneExpSampleFile = getRefDataFile(cmd, REF_RNA_GENE_EXP_SAMPLE_FILE, REF_FILE_GENE_EXP_SAMPLE);
        RefAltSjSampleFile = getRefDataFile(cmd, REF_RNA_ALT_SJ_SAMPLE_FILE, REF_FILE_ALT_SJ_SAMPLE);

        OutputDir = parseOutputDir(cmd);
    }

    private String getRefDataFile(final CommandLine cmd, final String configStr, final String defaultFilename)
    {
        final String fileName = cmd.getOptionValue(configStr, defaultFilename);
        return RefDataDir + fileName;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(REF_DATA_DIR, true, "Reference data directory");
        options.addOption(REF_SAMPLE_DATA_FILE, true, "Ref sample data file");
        options.addOption(REF_SNV_SAMPLE_POS_FREQ_FILE, true, "Ref SNV position frequency matrix data file");
        options.addOption(REF_SNV_COUNTS_FILE, true, "Ref SNV trinucleotide matrix data file");
        options.addOption(REF_RNA_GENE_EXP_SAMPLE_FILE, true, "Ref sample RNA gene expression cohort data file");
        options.addOption(REF_RNA_ALT_SJ_SAMPLE_FILE, true, "Ref sample RNA alternative SJ cohort data file");

        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
    }

}
