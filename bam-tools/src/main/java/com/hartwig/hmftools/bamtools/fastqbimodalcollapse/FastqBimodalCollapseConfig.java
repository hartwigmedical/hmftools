package com.hartwig.hmftools.bamtools.fastqbimodalcollapse;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE_DESC;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.common.bam.BamUtils.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.bamtools.tofastq.FileSplitMode;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.math3.exception.OutOfRangeException;
import org.jetbrains.annotations.Nullable;

public class FastqBimodalCollapseConfig
{
    public final String SampleId;
    public final String Fastq1Path;
    public final String Fastq2Path;
    public final String OutputDir;

    private static final String SAMPLE_ID = "sample";
    private static final String FASTQ1_PATH = "fastq1";
    private static final String FASTQ2_PATH = "fastq2";
    private static final String OUTPUT_DIR = "output_dir";

    public FastqBimodalCollapseConfig(final ConfigBuilder configBuilder)
    {
        SampleId = configBuilder.getValue(SAMPLE_ID);
        Fastq1Path = configBuilder.getValue(FASTQ1_PATH);
        Fastq2Path = configBuilder.getValue(FASTQ2_PATH);
        OutputDir = configBuilder.getValue(OUTPUT_DIR);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE_ID, true, "ID of the sample");
        configBuilder.addConfigItem(FASTQ1_PATH, true, "Path to the fastq file of first in pair reads");
        configBuilder.addConfigItem(FASTQ2_PATH, true, "Path to the fastq file of second in pair reads");
        configBuilder.addConfigItem(OUTPUT_DIR, true, "Directory for output files");
    }
}
