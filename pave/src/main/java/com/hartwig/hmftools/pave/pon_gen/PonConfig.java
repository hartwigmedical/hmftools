package com.hartwig.hmftools.pave.pon_gen;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificRegions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_FILE;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.pave.annotation.ClinvarAnnotation;

public class PonConfig
{
    public final String VcfPath;
    public final String OutputFilename;
    public final String ExistingPonFilename;

    public final int QualCutoff;
    public final int MqfCutoff;
    public final int MinSamples;
    public final boolean WriteFinal;
    public final boolean SkipGermlineIndelCheck;

    public final int RefSampleGenoptypeIndex;
    public final RefGenomeVersion RefGenVersion;
    public final int Threads;
    public final int PartitionSize;

    public final SpecificRegions SpecificChrRegions;

    private static final String VCF_PATH = "vcf_path";
    private static final String OUTPUT_PON_FILENAME = "output_pon_file";
    private static final String QUAL_CUTOFF = "qual_cutoff";
    private static final String MQF_CUTOFF = "mqf_cutoff";
    private static final String MIN_SAMPLES = "min_samples";
    private static final String REF_SAMPLE_GENOTYPE_INDEX = "ref_sample_index";
    private static final String PARTITION_SIZE = "partition_size";

    public static final String WRITE_FINAL_PON = "write_final";
    public static final String MANUAL_ENTRIES = "manual_entries";
    public static final String SOMATIC_HOTSPOT = "somatic_hotspots";
    public static final String GERMLINE_HOTSPOT = "germline_hotspots";
    public static final String SKIP_GERMLINE_INDEL_CHECK = "skip_germline_indel_check";

    // constants
    private static final int DEFAULT_MIN_SAMPLES = 3;
    protected static final int GERMLINE_CLINVAR_MIN_SAMPLES = 10;
    protected static final int GERMLINE_CLINVAR_MAX_REPEAT = 3;
    private static final int DEFAULT_MIN_MAP_QUAL = -10;
    private static final int DEFAULT_MIN_QUAL = 40;
    private static final int DEFAULT_PARTITION_SIZE = 10_000_000;

    public PonConfig(final ConfigBuilder configBuilder)
    {
        VcfPath = configBuilder.getValue(VCF_PATH);
        OutputFilename = configBuilder.getValue(OUTPUT_PON_FILENAME);

        QualCutoff = configBuilder.getInteger(QUAL_CUTOFF);
        MqfCutoff = configBuilder.getInteger(MQF_CUTOFF);
        MinSamples = configBuilder.getInteger(MIN_SAMPLES);
        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);
        RefSampleGenoptypeIndex = configBuilder.getInteger(REF_SAMPLE_GENOTYPE_INDEX);

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        WriteFinal = configBuilder.hasFlag(WRITE_FINAL_PON);
        SkipGermlineIndelCheck = configBuilder.hasFlag(SKIP_GERMLINE_INDEL_CHECK);

        ExistingPonFilename = configBuilder.getValue(PON_FILE);

        PV_LOGGER.info("key config: minSamples({}) cut-offs(qual={} mqf={}) skipGermlineIndelCheck({})",
                MinSamples, QualCutoff, MqfCutoff, SkipGermlineIndelCheck);

        Threads = parseThreads(configBuilder);

        SpecificChrRegions = SpecificRegions.from(configBuilder);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        addSampleIdFile(configBuilder, true);
        configBuilder.addConfigItem(VCF_PATH, true, "VCF path for samples");
        configBuilder.addConfigItem(OUTPUT_PON_FILENAME, false, "Output PON filename");
        configBuilder.addInteger(MIN_SAMPLES, "Min samples for variant to be included in PON", DEFAULT_MIN_SAMPLES);
        configBuilder.addInteger(QUAL_CUTOFF, "Qual cut-off for variant inclusion", DEFAULT_MIN_QUAL);
        configBuilder.addInteger(MQF_CUTOFF, "MQF cut-off for variant inclusion", DEFAULT_MIN_MAP_QUAL);
        configBuilder.addInteger(REF_SAMPLE_GENOTYPE_INDEX, "Ref sample genotype index", -1);
        configBuilder.addInteger(PARTITION_SIZE, "Partition size", DEFAULT_PARTITION_SIZE);
        configBuilder.addConfigItem(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        configBuilder.addConfigItem(MANUAL_ENTRIES, false, "Manual PON entries in form Chr:Pos:Ref:Alt separated by ';'");

        configBuilder.addFlag(WRITE_FINAL_PON, "Write final PON without annotations");
        configBuilder.addFlag(SKIP_GERMLINE_INDEL_CHECK, "Skip germline indel repeat filter");

        configBuilder.addPath(PON_FILE, false, "PON entries");
        ClinvarAnnotation.addConfig(configBuilder);
        addEnsemblDir(configBuilder, false);

        configBuilder.addPath(SOMATIC_HOTSPOT, false, "Somatic hotspot file");
        configBuilder.addPath(GERMLINE_HOTSPOT, false, "Germline hotspot file");

        addThreadOptions(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);

        addLoggingOptions(configBuilder);
    }
}
