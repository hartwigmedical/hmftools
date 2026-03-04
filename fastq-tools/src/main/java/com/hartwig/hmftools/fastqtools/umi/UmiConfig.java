package com.hartwig.hmftools.fastqtools.umi;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

public class UmiConfig
{
    public final String FastqFiles;
    public final String OutputDir;
    public final String OutputId;
    public final int UmiLength;
    public final String UmiDelim;
    public final int AdapterLength;
    public final int AdapterUmiLength;
    public final String AdapterSequence;

    public final String KnownUmiFile;
    public final int KnownUmiBaseDiff;

    public final String AdapterSequenceReversed;

    private static final String FASTQ_FILES = "fastq_files";
    private static final String UMI_LENGTH = "umi_length";
    private static final String UMI_DELIM = "umi_delim";
    private static final String ADAPTER_LENGTH = "adapter_length";
    private static final String ADAPTER_SEQUENCE = "adapter_seq";
    private static final String KNOWN_UMI_FILE = "known_umi_file";
    private static final String KNOWN_UMI_BASE_DIFF = "known_umi_base_diff";

    private static final int DEFAULT_KNOWN_UMI_BASE_DIFF = 1;

    public UmiConfig(final ConfigBuilder configBuilder)
    {
        FastqFiles = configBuilder.getValue(FASTQ_FILES);
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID, "umi");

        UmiLength = configBuilder.getInteger(UMI_LENGTH);

        AdapterLength = configBuilder.getInteger(ADAPTER_LENGTH);
        AdapterSequence = configBuilder.getValue(ADAPTER_SEQUENCE);

        KnownUmiFile = configBuilder.getValue(KNOWN_UMI_FILE);
        KnownUmiBaseDiff = configBuilder.getInteger(KNOWN_UMI_BASE_DIFF);

        UmiDelim = configBuilder.getValue(UMI_DELIM);

        if(AdapterSequence != null)
        {
            AdapterUmiLength = AdapterSequence != null ? AdapterSequence.length() + UmiLength : 0;
            AdapterSequenceReversed = AdapterSequence != null ? Nucleotides.reverseComplementBases(AdapterSequence) : null;
        }
        else
        {
            AdapterUmiLength = AdapterLength > 0 ? UmiLength + AdapterLength : UmiLength;
            AdapterSequenceReversed = null;
        }
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(FASTQ_FILES, true, "Fastq file-pair path, separated by delim ';'");
        configBuilder.addPath(KNOWN_UMI_FILE, false, "File with known UMI sequences");
        configBuilder.addInteger(UMI_LENGTH,"UMI length", 0);
        configBuilder.addConfigItem(ADAPTER_SEQUENCE, "Adapter sequence (optional)");
        configBuilder.addConfigItem(UMI_DELIM, "UMI delimiter");
        configBuilder.addInteger(KNOWN_UMI_BASE_DIFF, "Max permitted base difference for known UMI", DEFAULT_KNOWN_UMI_BASE_DIFF);
        configBuilder.addInteger(ADAPTER_LENGTH, "Adapter length", 0);

        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
