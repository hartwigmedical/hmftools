package com.hartwig.hmftools.svtools.pon;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.svtools.pon.PonBuilder.PON_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import com.google.common.collect.Lists;

public class PonConfig
{
    public final List<String> SampleIds;
    public final String OutputDir;
    public final List<String> VcfFilePatterns;
    public final int MinPonWriteCount;
    public final int Threads;

    final Map<String,String> SampleVcfFiles;

    private static final String SAMPLE_VCFS_FILE = "sample_vcfs_file";
    private static final String VCF_FILE_PATTERNS = "vcf_file_patterns";
    private static final String MIN_PON_WRITE_COUNT = "min_pon_write_count";

    public PonConfig(final ConfigBuilder configBuilder)
    {
        SampleIds = Lists.newArrayList();
        OutputDir = parseOutputDir(configBuilder);
        MinPonWriteCount = configBuilder.getInteger(MIN_PON_WRITE_COUNT);
        Threads = parseThreads(configBuilder);

        SampleVcfFiles = Maps.newHashMap();

        if(configBuilder.hasFlag(SAMPLE_VCFS_FILE))
        {
            VcfFilePatterns = Lists.newArrayList();
            populateSampleVcfFilepaths(configBuilder.getValue(SAMPLE_VCFS_FILE));
            SampleVcfFiles.keySet().forEach(x -> SampleIds.add(x));
        }
        else
        {
            SampleIds.addAll(ConfigUtils.loadSampleIdsFile(configBuilder));
            VcfFilePatterns = Arrays.stream(configBuilder.getValue(VCF_FILE_PATTERNS).split(";")).collect(Collectors.toList());
        }
    }

    private void populateSampleVcfFilepaths(final String filename)
    {
        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
                return;

            String header = fileContents.get(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");
            fileContents.remove(0);
            int sampleIdIndex = fieldsIndexMap.get("SampleId");
            int vcfIndex = fieldsIndexMap.get("VcfFile");

            for(String line : fileContents)
            {
                String[] values = line.split(",", -1);
                SampleVcfFiles.put(values[sampleIdIndex], values[vcfIndex]);
            }

            PON_LOGGER.info("loaded {} sample-VCF file entries", SampleVcfFiles.size());
        }
        catch (IOException e)
        {
            PON_LOGGER.error("failed to read file({}): {}", filename, e.toString());
        }
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        addSampleIdFile(configBuilder, true);
        configBuilder.addConfigItem(SAMPLE_VCFS_FILE, true, "CSV file with 'SampleId,VcfFile' locations");
        configBuilder.addConfigItem(VCF_FILE_PATTERNS, true, "VCF file IDs, eg 'gridss.vcf' separated by ';'");
        configBuilder.addInteger(MIN_PON_WRITE_COUNT, "Min observations of SV or SGL to include in PON", 2);
        addThreadOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
        addOutputDir(configBuilder);
    }
}
