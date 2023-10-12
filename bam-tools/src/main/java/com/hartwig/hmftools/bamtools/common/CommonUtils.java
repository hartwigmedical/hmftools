package com.hartwig.hmftools.bamtools.common;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificChromsomesOrRegions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.bamtools.metrics.MetricsConfig;
import com.hartwig.hmftools.common.genome.bed.BedFileReader;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class CommonUtils
{
    public static final String APP_NAME = "BamTools";

    // constants
    public static final int DEFAULT_CHR_PARTITION_SIZE = 1000000;
    public static final int DEFAULT_READ_LENGTH = 151;

    // config strings
    public static final String BAM_FILE = "bam_file";
    public static final String PARTITION_SIZE = "partition_size";
    public static final String REGIONS_FILE = "regions_file";
    public static final String READ_LENGTH = "read_length";

    public static final String BAM_FILE_TYPE = "bam";
    public static final String FLAGSTAT_FILE_TYPE = "flagstat";

    public static final Logger BT_LOGGER = LogManager.getLogger(MetricsConfig.class);

    public static void addCommonCommandOptions(final ConfigBuilder configBuilder)
    {
        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
        addRefGenomeConfig(configBuilder, true);;
        addSpecificChromosomesRegionsConfig(configBuilder);

        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        configBuilder.addPath(BAM_FILE, true, "BAM file location");

        configBuilder.addPath(
                REGIONS_FILE, false,
                "TSV or BED file with regions to analyse, expected columns Chromosome,PositionStart,PositionEnd or no headers");

        configBuilder.addInteger(READ_LENGTH, "Read length", DEFAULT_READ_LENGTH);
    }

    public static boolean loadSpecificRegionsConfig(
            final ConfigBuilder configBuilder, final List<String> specificChromosomes, final List<ChrBaseRegion> specificRegions)
    {
        if(configBuilder.isRegistered(REGIONS_FILE) && configBuilder.hasValue(REGIONS_FILE))
        {
            return BedFileReader.loadBedFile(configBuilder.getValue(REGIONS_FILE), specificRegions);
        }
        else
        {
            try
            {
                loadSpecificChromsomesOrRegions(configBuilder, specificChromosomes, specificRegions);
            }
            catch(Exception e)
            {
                BT_LOGGER.error("failed to load specific regions: {}", e.toString());
                return false;
            }
        }

        return true;
    }

    public static boolean checkFileExists(final String filename)
    {
        if(!Files.exists(Paths.get(filename)))
        {
            BT_LOGGER.error("invalid file path: {}", filename);
            return false;
        }

        return true;
    }

    public static String formFilename(
            final String sampleId, final String bamFile, final String outputDir, final String outputId, final String fileType)
    {
        String filename;

        if(sampleId != null && !sampleId.isEmpty())
        {
            filename = outputDir + sampleId;
        }
        else
        {
            filename = bamFile.substring(0, bamFile.indexOf(".bam"));
        }

        if(fileType.equals(FLAGSTAT_FILE_TYPE))
        {
            filename += ".flagstat";
        }
        else if(!fileType.equals(BAM_FILE_TYPE))
        {
            filename += ".bam_" + fileType;
        }

        if(outputId != null)
        {
            filename += "." + outputId;
        }

        if(fileType.equals(BAM_FILE_TYPE))
        {
            filename += ".bam";
        }
        else if(!fileType.equals(FLAGSTAT_FILE_TYPE))
        {
            filename += TSV_EXTENSION;
        }

        return filename;
    }
}
