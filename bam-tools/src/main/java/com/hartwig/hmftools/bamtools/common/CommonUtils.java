package com.hartwig.hmftools.bamtools.common;

import static com.hartwig.hmftools.bamtools.metrics.MetricsConfig.BT_LOGGER;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificChromsomesOrRegions;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.common.genome.bed.BedFileReader;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public final class CommonUtils
{
    // constants
    public static final int DEFAULT_CHR_PARTITION_SIZE = 1000000;

    // config strings
    public static final String SAMPLE = "sample";
    public static final String BAM_FILE = "bam_file";
    public static final String PARTITION_SIZE = "partition_size";
    public static final String REGIONS_BED_FILE = "regions_bed_file";

    public static void addCommonCommandOptions(final Options options)
    {
        addOutputOptions(options);
        addLoggingOptions(options);
        addThreadOptions(options);
        addRefGenomeConfig(options);;
        addSpecificChromosomesRegionsConfig(options);

        options.addOption(SAMPLE, true, "Tumor sample ID");
        options.addOption(BAM_FILE, true, "BAM file location");
        options.addOption(REGIONS_BED_FILE, true, "BED file with regions to analyse");
    }

    public static boolean loadSpecificRegionsConfig(
            final CommandLine cmd, final List<String> specificChromosomes, final List<ChrBaseRegion> specificRegions)
    {
        if(cmd.hasOption(REGIONS_BED_FILE))
        {
            return BedFileReader.loadBedFile(cmd.getOptionValue(REGIONS_BED_FILE), specificRegions);
        }
        else
        {
            try
            {
                loadSpecificChromsomesOrRegions(cmd, specificChromosomes, specificRegions, BT_LOGGER);
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

    public static String formFilename(final String sampleId, final String outputDir, final String outputId, final String fileType)
    {
        String filename = outputDir + sampleId;

        filename += ".bam_" + fileType;

        if(outputId != null)
            filename += "." + outputId;

        filename += ".csv";

        return filename;
    }
}
