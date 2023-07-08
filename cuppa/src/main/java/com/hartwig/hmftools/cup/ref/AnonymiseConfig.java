package com.hartwig.hmftools.cup.ref;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
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

import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

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

    public AnonymiseConfig(final ConfigBuilder configBuilder)
    {
        RefDataDir = checkAddDirSeparator(configBuilder.getValue(REF_DATA_DIR, ""));

        RefSampleDataFile = getRefDataFile(configBuilder, REF_SAMPLE_DATA_FILE, REF_FILE_SAMPLE_DATA);
        RefSnvCountsFile = getRefDataFile(configBuilder, REF_SNV_COUNTS_FILE, REF_FILE_SNV_COUNTS);
        RefSnvSamplePosFreqFile = getRefDataFile(configBuilder, REF_SNV_SAMPLE_POS_FREQ_FILE, REF_FILE_SAMPLE_POS_FREQ_COUNTS);
        RefGeneExpSampleFile = getRefDataFile(configBuilder, REF_RNA_GENE_EXP_SAMPLE_FILE, REF_FILE_GENE_EXP_SAMPLE, true);
        RefAltSjSampleFile = getRefDataFile(configBuilder, REF_RNA_ALT_SJ_SAMPLE_FILE, REF_FILE_ALT_SJ_SAMPLE, true);

        OutputDir = parseOutputDir(configBuilder);
    }

    private String getRefDataFile(final ConfigBuilder configBuilder, final String configStr, final String defaultFilename)
    {
        return getRefDataFile(configBuilder, configStr, defaultFilename, false);
    }

    private String getRefDataFile(final ConfigBuilder configBuilder, final String configStr, final String defaultFilename, boolean checkZipped)
    {
        final String fileName = RefDataDir + configBuilder.getValue(configStr, defaultFilename);

        if(checkZipped && !Files.exists(Paths.get(fileName)) && Files.exists(Paths.get(fileName + ".gz")))
            return fileName + ".gz";

        return fileName;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(REF_DATA_DIR, true, "Reference data directory");
        configBuilder.addPath(REF_SAMPLE_DATA_FILE, true, "Ref sample data file");
        configBuilder.addPath(REF_SNV_SAMPLE_POS_FREQ_FILE, false, "Ref SNV position frequency matrix data file");
        configBuilder.addPath(REF_SNV_COUNTS_FILE, false, "Ref SNV trinucleotide matrix data file");
        configBuilder.addPath(REF_RNA_GENE_EXP_SAMPLE_FILE, false, "Ref sample RNA gene expression cohort data file");
        configBuilder.addPath(REF_RNA_ALT_SJ_SAMPLE_FILE, false, "Ref sample RNA alternative SJ cohort data file");

        addLoggingOptions(configBuilder);
        addOutputOptions(configBuilder);
    }
}
