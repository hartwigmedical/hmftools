package com.hartwig.hmftools.wisp.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.COBALT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.COBALT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.SampleData.ctDnaSamplesFromStr;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.DEFAULT_NOISE_READS_PER_MILLION;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.DEFAULT_NOISE_READS_PER_MILLION_DUAL_STRAND;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.wisp.purity.variant.ProbeVariantCache;

public class PurityConfig
{
    public final List<SampleData> Samples;

    public final List<PurityMethod> PurityMethods;

    public final String SomaticVcf;
    public final String SomaticDir; // if different from sample data dir and not specifying a VCF path
    public final String SampleDataDir;
    public final String PurpleDir;
    public final String CobaltDir;

    public final ProbeVariantCache ProbeVariants;

    public final String OutputDir;
    public final String PlotDir;
    public final String OutputId;
    public final RefGenomeInterface RefGenome;
    public final Set<WriteType> WriteTypes;
    public final double NoiseReadsPerMillion;
    public final double NoiseReadsPerMillionDualStrand;
    public final double GcRatioMin;
    public final int Threads;

    private static final String PATIENT_ID = "patient_id";
    private static final String TUMOR_ID = "tumor_id";
    private static final String CTDNA_SAMPLES = "ctdna_samples";
    private static final String PURITY_METHODS = "purity_methods";
    private static final String SOMATIC_VCF = "somatic_vcf";
    private static final String SOMATIC_DIR = "somatic_dir";
    private static final String PLOT_DIR = "plot_dir";
    private static final String NOISE_READS_PER_MILLION = "noise_per_mill";
    private static final String NOISE_READS_PER_MILLION_DUAL = "noise_per_mill_dual";
    private static final String GC_RATIO_MIN = "gc_ratio_min";
    private static final String WRITE_TYPES = "write_types";
    private static final String PROBE_VARIANTS_FILE = "probe_variants_file";

    public PurityConfig(final ConfigBuilder configBuilder)
    {
        SampleDataDir = checkAddDirSeparator(configBuilder.getValue(SAMPLE_DATA_DIR_CFG));

        GcRatioMin = configBuilder.getDecimal(GC_RATIO_MIN);

        Samples = Lists.newArrayList();
        loadSampleData(configBuilder);

        PurityMethods = Lists.newArrayList();

        if(configBuilder.hasValue(PURITY_METHODS))
        {
            Arrays.stream(configBuilder.getValue(PURITY_METHODS).split(ITEM_DELIM, -1))
                    .forEach(x -> PurityMethods.add(PurityMethod.valueOf(x)));
        }
        else
        {
            Arrays.stream(PurityMethod.values()).forEach(x -> PurityMethods.add(x));
        }

        SomaticVcf = configBuilder.getValue(SOMATIC_VCF);
        SomaticDir = checkAddDirSeparator(configBuilder.getValue(SOMATIC_DIR, SampleDataDir));
        PurpleDir = checkAddDirSeparator(configBuilder.getValue(PURPLE_DIR_CFG, SampleDataDir));
        CobaltDir = checkAddDirSeparator(configBuilder.getValue(COBALT_DIR_CFG, SampleDataDir));
        OutputDir = checkAddDirSeparator(configBuilder.getValue(OUTPUT_DIR, SampleDataDir));
        OutputId = configBuilder.getValue(OUTPUT_ID);

        PlotDir = checkAddDirSeparator(configBuilder.getValue(PLOT_DIR, OutputDir));

        ProbeVariants = new ProbeVariantCache(configBuilder.getValue(PROBE_VARIANTS_FILE));

        RefGenome = configBuilder.hasValue(REF_GENOME) ? loadRefGenome(configBuilder.getValue(REF_GENOME)) : null;

        NoiseReadsPerMillion = configBuilder.getDecimal(NOISE_READS_PER_MILLION);
        NoiseReadsPerMillionDualStrand = configBuilder.getDecimal(NOISE_READS_PER_MILLION_DUAL);

        WriteTypes = Sets.newHashSet();

        if(configBuilder.hasValue(WRITE_TYPES))
        {
            String writeTypes = configBuilder.getValue(WRITE_TYPES);

            if(writeTypes.equals(WriteType.ALL))
                Arrays.stream(WriteType.values()).forEach(x -> WriteTypes.add(x));
            else
                Arrays.stream(writeTypes.split(ITEM_DELIM, -1)).forEach(x -> WriteTypes.add(WriteType.valueOf(x)));
        }

        Threads = parseThreads(configBuilder);
    }

    public boolean writeType(final WriteType writeType) { return WriteTypes.contains(writeType); }
    public boolean hasSyntheticTumor() { return PurpleDir == null || PurpleDir.isEmpty(); }
    public boolean multiplePatients() { return Samples.size() > 1; }
    public boolean multipleSamples() { return multiplePatients() || Samples.stream().mapToInt(x -> x.CtDnaSamples.size()).sum() > 1; }
    public boolean hasBatchControls() { return Samples.stream().anyMatch(x -> x.isBatchControl()); }

    public String getPurpleDir(final String sampleId) { return convertWildcardSamplePath(PurpleDir, sampleId); }
    public String getSomaticVcf(final String sampleId) { return convertWildcardSamplePath(SomaticVcf, sampleId); }
    public String getCobaltDir(final String sampleId) { return convertWildcardSamplePath(CobaltDir, sampleId); }

    public double noiseRate(boolean useDual)
    {
        return (useDual ? NoiseReadsPerMillionDualStrand : NoiseReadsPerMillion) / 1_000_000d;
    }

    private void loadSampleData(final ConfigBuilder configBuilder)
    {
        if(configBuilder.hasValue(SAMPLE_ID_FILE))
        {
            String filename = configBuilder.getValue(SAMPLE_ID_FILE);

            if(!Files.exists(Paths.get(filename)))
                filename = SampleDataDir + filename;

            Samples.addAll(SampleData.loadSampleDataFile(filename));
        }
        else
        {
            Samples.add(new SampleData(
                    configBuilder.getValue(PATIENT_ID),
                    configBuilder.getValue(TUMOR_ID),
                    ctDnaSamplesFromStr(configBuilder.getValue(CTDNA_SAMPLES)), "", GcRatioMin > 0));
        }

        CT_LOGGER.info("loaded {} patients and {} ctDNA samples",
                Samples.size(), Samples.stream().mapToInt(x -> x.CtDnaSamples.size()).sum());
    }

    public String formFilename(final String fileType)
    {
        String fileName = OutputDir;

        if(multiplePatients())
        {
            fileName += "wisp_cohort.";
        }
        else if(Samples.get(0).CtDnaSamples.size() > 1)
        {
            fileName += format("%s.wisp.", Samples.get(0).PatientId);
        }
        else
        {
            fileName += format("%s_%s.wisp.", Samples.get(0).PatientId, Samples.get(0).CtDnaSamples.get(0));
        }

        fileName += fileType;

        if(OutputId != null)
            fileName += "." + OutputId;

        fileName += TSV_EXTENSION;

        return fileName;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE_ID_FILE, false, "Patient and sample data file: PatientId,TumorId,CtDnaSampleIds");
        configBuilder.addConfigItem(PATIENT_ID, false, "Patient ID");
        configBuilder.addConfigItem(TUMOR_ID, false, "Original tumor sample ID");
        configBuilder.addConfigItem(CTDNA_SAMPLES, false, "List of ctDNA sample IDs separated by ','");

        configBuilder.addConfigItem(
                PURITY_METHODS, false,
                "List of purity methods separated by ',' default(all) from: "
                        + Arrays.stream(PurityMethod.values()).map(x -> x.toString()).collect(Collectors.joining(",")));

        configBuilder.addConfigItem(SOMATIC_VCF, false, "Somatic VCF files, separated by ','", "");
        configBuilder.addConfigItem(SOMATIC_DIR, false, "Somatic VCF directory");
        configBuilder.addConfigItem(SAMPLE_DATA_DIR_CFG, false, SAMPLE_DATA_DIR_DESC);
        configBuilder.addConfigItem(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addConfigItem(COBALT_DIR_CFG, false, COBALT_DIR_DESC);
        configBuilder.addConfigItem(PLOT_DIR, false, "Plot output directory, defaults to sample or output dir");

        configBuilder.addConfigItem(
                WRITE_TYPES, "Output file types: default(none), 'ALL' or set separated by ',': "
                        + Arrays.stream(WriteType.values()).map(x -> x.toString()).collect(Collectors.joining(",")));

        configBuilder.addPath(REF_GENOME, true, REF_GENOME_CFG_DESC);
        configBuilder.addPath(PROBE_VARIANTS_FILE, false, "File defining the probe variants");

        configBuilder.addDecimal(
                NOISE_READS_PER_MILLION, "Expected reads-per-million from noise", DEFAULT_NOISE_READS_PER_MILLION);

        configBuilder.addDecimal(
                NOISE_READS_PER_MILLION_DUAL,
                "Expected reads-per-million from noise for dual-strand reads", DEFAULT_NOISE_READS_PER_MILLION_DUAL_STRAND);

        configBuilder.addRequiredDecimal(GC_RATIO_MIN,"GC ratio minimum permitted");

        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
    }
}
