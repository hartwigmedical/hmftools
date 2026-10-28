package com.hartwig.hmftools.sigs.fitter;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SAMPLE_COUNTS_FILE;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SIGNATURES_FILE;
import static com.hartwig.hmftools.sigs.loaders.PositionFreqBuilder.POSITION_BUCKET_SIZE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class SigsConfig
{
    public final RefGenomeVersion RefGenVersion;
    public final List<String> SampleIdList;
    public final String SnvCountsFile;
    public final String SignaturesFile;
    public final String VcfFile;
    public final String OutputDir;
    public final String OutputId;

    public final boolean WritePosFreqCoords;
    public final int PositionBucketSize;
    public final double MinFit;
    public final double MinFitPerc;

    private static final double DEFAULT_MIN_ALLOCATION = 0.01;
    private static final double DEFAULT_MIN_ALLOCATION_PERC = 0.0005;

    public static final String SOMATIC_VCF_FILE = "somatic_vcf_file";
    public static final String MIN_ALLOC = "min_alloc";
    public static final String MIN_ALLOC_PERC = "min_alloc_perc";
    public static final String WRITE_POS_COORDS = "write_pos_freq_coords";

    public SigsConfig(final ConfigBuilder configBuilder)
    {
        SnvCountsFile = configBuilder.getValue(SAMPLE_COUNTS_FILE);
        SignaturesFile = configBuilder.getValue(SIGNATURES_FILE);
        VcfFile = configBuilder.getValue(SOMATIC_VCF_FILE);

        SampleIdList = Lists.newArrayList();

        if(configBuilder.hasValue(SAMPLE))
        {
            SampleIdList.add(configBuilder.getValue(SAMPLE));
        }
        else if(configBuilder.hasValue(SAMPLE_ID_FILE))
        {
            SampleIdList.addAll(loadSampleIdsFile(configBuilder));
        }

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        WritePosFreqCoords = configBuilder.hasFlag(WRITE_POS_COORDS);

        PositionBucketSize = configBuilder.getInteger(POSITION_BUCKET_SIZE);
        MinFit = configBuilder.getDecimal(MIN_ALLOC);
        MinFitPerc = configBuilder.getDecimal(MIN_ALLOC_PERC);

    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        configBuilder.addPath(SAMPLE_ID_FILE, false, SAMPLE_ID_FILE_DESC);

        // either a VCF for a single sample or a sample matrix file
        configBuilder.addPath(SOMATIC_VCF_FILE, false, "Somatic variant VCF file");
        configBuilder.addPath(SAMPLE_COUNTS_FILE, false, "Path to the main input file");

        configBuilder.addPath(SIGNATURES_FILE, false, "Signature definitions");
        addRefGenomeVersion(configBuilder);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.addInteger(POSITION_BUCKET_SIZE, "Position bucket size", 0);

        configBuilder.addDecimal(MIN_ALLOC_PERC, "Min signature allocation as percentage", DEFAULT_MIN_ALLOCATION_PERC);
        configBuilder.addDecimal(MIN_ALLOC, "Min signature allocation", DEFAULT_MIN_ALLOCATION);
        configBuilder.addFlag(WRITE_POS_COORDS, "Include coordinates with pos-frequency counts");
    }
}
