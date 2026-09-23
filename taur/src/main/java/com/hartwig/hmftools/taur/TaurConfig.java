package com.hartwig.hmftools.taur;

import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.taur.FastaCommon.FASTQ_SUFFIX_SHORT;
import static com.hartwig.hmftools.taur.FastaCommon.FASTQ_SUFFIX_STANDARD;
import static com.hartwig.hmftools.taur.FastaCommon.FASTQ_ZIP_EXTENSION;

import java.io.File;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

public class TaurConfig
{
    public final String FastqFiles;
    public final String OutputDir;
    public final String OutputId;
    public final int Threads;
    public final int ReadsPerChunk;

    public final int UmiLength;
    public final String UmiDelim;
    public final int AdapterLength;
    public final String AdapterSequence;

    public final String KnownUmiFile;
    public final int KnownUmiBaseDiff;
    public final boolean KnownUmiUseNumeric;

    private static final String FASTQ_FILES = "fastq_files";
    private static final String UMI_LENGTH = "umi_length";
    private static final String UMI_DELIM = "umi_delim";
    private static final String ADAPTER_LENGTH = "adapter_length";
    private static final String ADAPTER_SEQUENCE = "adapter_seq";
    private static final String CHUNK_SIZE = "chunk_size";

    private static final String KNOWN_UMI_FILE = "known_umi_file";
    private static final String KNOWN_UMI_BASE_DIFF = "known_umi_base_diff";
    private static final String KNOWN_UMI_USE_NUMERIC = "known_umi_use_numeric";

    protected static final String FASTQ_FILES_DELIM = ";";

    private static final int DEFAULT_KNOWN_UMI_BASE_DIFF = 0;
    private static final int READ_GROUPS_PER_CHUNK_SIZE_IN_BYTES = 100_000;

    public TaurConfig(final ConfigBuilder configBuilder)
    {
        FastqFiles = configBuilder.getValue(FASTQ_FILES);
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID, "umi");

        Threads = parseThreads(configBuilder);
        ReadsPerChunk = configBuilder.getInteger(CHUNK_SIZE);

        UmiLength = configBuilder.getInteger(UMI_LENGTH);

        AdapterLength = configBuilder.getInteger(ADAPTER_LENGTH);
        AdapterSequence = configBuilder.getValue(ADAPTER_SEQUENCE);

        KnownUmiFile = configBuilder.getValue(KNOWN_UMI_FILE);
        KnownUmiBaseDiff = configBuilder.getInteger(KNOWN_UMI_BASE_DIFF);
        KnownUmiUseNumeric = configBuilder.hasFlag(KNOWN_UMI_USE_NUMERIC);

        UmiDelim = configBuilder.getValue(UMI_DELIM);
   }

   public String outputFastqFilename(final String fastqFilename)
   {
       String fastqFile = fastqFilename.substring(fastqFilename.lastIndexOf(File.separator) + 1);

       int extensionIndex = fastqFile.contains(FASTQ_SUFFIX_STANDARD) ?
               fastqFile.lastIndexOf("." + FASTQ_SUFFIX_STANDARD) : fastqFile.lastIndexOf("." + FASTQ_SUFFIX_SHORT);

       return OutputDir + fastqFile.substring(0, extensionIndex) + '.' + OutputId + FASTQ_ZIP_EXTENSION;
   }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(FASTQ_FILES, true, "Fastq file-pair path, separated by delim ';'");
        configBuilder.addInteger(UMI_LENGTH,"UMI length", 0);
        configBuilder.addConfigItem(ADAPTER_SEQUENCE, "Adapter sequence (optional)");
        configBuilder.addConfigItem(UMI_DELIM, "UMI delimiter");

        configBuilder.addPath(KNOWN_UMI_FILE, false, "File with known UMI sequences");
        configBuilder.addInteger(KNOWN_UMI_BASE_DIFF, "Max permitted base difference for known UMI", DEFAULT_KNOWN_UMI_BASE_DIFF);
        configBuilder.addFlag(KNOWN_UMI_USE_NUMERIC, "Translate known UMI to numerics, eg ABC+WXYZ -> 3+5");

        configBuilder.addInteger(ADAPTER_LENGTH, "Adapter length", 0);

        configBuilder.addInteger(CHUNK_SIZE,"Reads per size", READ_GROUPS_PER_CHUNK_SIZE_IN_BYTES);
        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
