package com.hartwig.hmftools.svtools.pon;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.svtools.pon.PonBuilder.PON_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.ConfigUtils;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;

public class PonConfig
{
    public final List<String> SampleIds;
    public final String OutputDir;
    public final List<String> VcfFilePatterns;
    public final int MinPonWriteCount;
    public final int Threads;

    final Map<String,String> SampleVcfFiles;

    private static final String SAMPLES_FILE = "samples_file";
    private static final String SAMPLE_VCFS_FILE = "sample_vcfs_file";
    private static final String VCF_FILE_PATTERNS = "vcf_file_patterns";
    private static final String MIN_PON_WRITE_COUNT = "min_pon_write_count";
    private static final String THREADS = "threads";

    public PonConfig(final CommandLine cmd)
    {
        SampleIds = Lists.newArrayList();
        OutputDir = parseOutputDir(cmd);
        MinPonWriteCount = Integer.parseInt(cmd.getOptionValue(MIN_PON_WRITE_COUNT, "2"));
        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "0"));

        SampleVcfFiles = Maps.newHashMap();

        if(cmd.hasOption(SAMPLE_VCFS_FILE))
        {
            VcfFilePatterns = Lists.newArrayList();
            populateSampleVcfFilepaths(cmd.getOptionValue(SAMPLE_VCFS_FILE));
            SampleVcfFiles.keySet().forEach(x -> SampleIds.add(x));
        }
        else
        {
            SampleIds.addAll(ConfigUtils.loadSampleIdsFile(cmd.getOptionValue(SAMPLES_FILE)));
            VcfFilePatterns = Arrays.stream(cmd.getOptionValue(VCF_FILE_PATTERNS).split(";")).collect(Collectors.toList());
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

    public static void addOptions(final Options options)
    {
        options.addOption(SAMPLES_FILE, true, "SampleIds file");
        options.addOption(SAMPLE_VCFS_FILE, true, "CSV file with 'SampleId,VcfFile' locations");
        options.addOption(VCF_FILE_PATTERNS, true, "VCF file IDs, eg 'gridss.vcf' separated by ';'");
        options.addOption(MIN_PON_WRITE_COUNT, true, "Min observations of SV or SGL to include in PON");
        options.addOption(THREADS, true, "Thread count (default: 0, none)");
        ConfigUtils.addLoggingOptions(options);
        addOutputDir(options);
    }
}
