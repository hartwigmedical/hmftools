package com.hartwig.hmftools.svtools.pon;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.ConfigUtils;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class PonConfig
{
    public final List<String> SampleIds;
    public final String RootDirectory;
    public final String OutputDir;
    public final List<String> VcfFileIds;
    public final int MinPonWriteCount;

    public static final String SAMPLES_FILE = "samples_file";
    public static final String ROOT_DIR = "root_dir";
    public static final String VCF_FILE_IDS = "vcf_file_ids";
    public static final String MIN_PON_WRITE_COUNT = "min_pon_write_count";

    public PonConfig(final CommandLine cmd)
    {
        SampleIds = ConfigUtils.loadSampleIdsFile(cmd.getOptionValue(SAMPLES_FILE));
        RootDirectory = checkAddDirSeparator(cmd.getOptionValue(ROOT_DIR));
        OutputDir = parseOutputDir(cmd);
        MinPonWriteCount = Integer.parseInt(MIN_PON_WRITE_COUNT);

        VcfFileIds = Arrays.stream(cmd.getOptionValue(VCF_FILE_IDS).split(";")).collect(Collectors.toList());
    }

    public static void addOptions(final Options options)
    {
        options.addOption(SAMPLES_FILE, true, "SampleIds file");
        options.addOption(ROOT_DIR, true, "Root directory for sample files");
        options.addOption(VCF_FILE_IDS, true, "VCF file IDs, eg 'gridss.vcf' separated by ';'");
        options.addOption(MIN_PON_WRITE_COUNT, true, "Min observations of SV or SGL to include in PON");
        ConfigUtils.addLoggingOptions(options);
        addOutputDir(options);
    }
}
