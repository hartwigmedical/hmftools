package com.hartwig.hmftools.ctdna.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.COBALT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.COBALT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.DEFAULT_GC_RATIO_MIN;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.DEFAULT_NOISE_READS_PER_MILLION;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.DEFAULT_NOISE_READS_PER_MILLION_DUAL_STRAND;
import static com.hartwig.hmftools.ctdna.common.SampleData.ctDnaSamplesFromStr;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.ctdna.common.SampleData;

public class PurityConfig
{
    public final List<SampleData> Samples;

    public final List<PurityMethod> PurityMethods;

    public final String SomaticVcf;
    public final String SampleDataDir;
    public final String PurpleDir;
    public final String CobaltDir;
    public final String OutputDir;
    public final String OutputId;
    public final RefGenomeInterface RefGenome;
    public final boolean WriteSomatics;
    public final boolean WriteCnRatios;
    public final boolean PlotCnFit;
    public final boolean ApplyDropout;
    public final boolean WriteFilteredSomatics;
    public final double NoiseReadsPerMillion;
    public final double NoiseReadsPerMillionDualStrand;
    public final double GcRatioMin;
    public final int Threads;

    private static final String PATIENT_ID = "patient_id";
    private static final String TUMOR_ID = "tumor_id";
    private static final String CTDNA_SAMPLES = "ctdna_samples";
    private static final String PURITY_METHODS = "purity_methods";
    private static final String SOMATIC_VCF = "somatic_vcf";
    private static final String WRITE_VARIANTS = "write_somatics";
    private static final String INCLUDE_FILTERED_VARIANTS = "write_filtered_somatics";
    private static final String WRITE_CN_RATIOS = "write_cn_ratios";
    private static final String PLOT_CN = "plot_cn_fit";
    private static final String NOISE_READS_PER_MILLION = "noise_per_mill";
    private static final String NOISE_READS_PER_MILLION_DUAL = "noise_per_mill_dual";
    private static final String APPLY_DROPOUT = "apply_dropout";
    private static final String GC_RATIO_MIN = "gc_ratio_min";

    public PurityConfig(final ConfigBuilder configBuilder)
    {
        SampleDataDir = checkAddDirSeparator(configBuilder.getValue(SAMPLE_DATA_DIR_CFG));

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
        PurpleDir = checkAddDirSeparator(configBuilder.getValue(PURPLE_DIR_CFG, SampleDataDir));
        CobaltDir = checkAddDirSeparator(configBuilder.getValue(COBALT_DIR_CFG, SampleDataDir));
        OutputDir = checkAddDirSeparator(configBuilder.getValue(OUTPUT_DIR, SampleDataDir));
        OutputId = configBuilder.getValue(OUTPUT_ID);

        RefGenome = configBuilder.hasValue(REF_GENOME) ? loadRefGenome(configBuilder.getValue(REF_GENOME)) : null;

        NoiseReadsPerMillion = configBuilder.getDecimal(NOISE_READS_PER_MILLION);
        NoiseReadsPerMillionDualStrand = configBuilder.getDecimal(NOISE_READS_PER_MILLION_DUAL);
        GcRatioMin = configBuilder.getDecimal(GC_RATIO_MIN);

        ApplyDropout = configBuilder.hasFlag(APPLY_DROPOUT);
        WriteSomatics = configBuilder.hasFlag(WRITE_VARIANTS);
        WriteCnRatios = configBuilder.hasFlag(WRITE_CN_RATIOS);
        WriteFilteredSomatics = configBuilder.hasFlag(INCLUDE_FILTERED_VARIANTS);
        PlotCnFit = configBuilder.hasFlag(PLOT_CN);
        Threads = parseThreads(configBuilder);
    }

    private void loadSampleData(final ConfigBuilder configBuilder)
    {
        if(configBuilder.hasValue(SAMPLE_ID_FILE))
        {
            Samples.addAll(SampleData.loadSampleDataFile(configBuilder.getValue(SAMPLE_ID_FILE)));
        }
        else
        {
            Samples.add(new SampleData(
                    configBuilder.getValue(PATIENT_ID),
                    configBuilder.getValue(TUMOR_ID),
                    ctDnaSamplesFromStr(configBuilder.getValue(CTDNA_SAMPLES)), ""));
        }

        CT_LOGGER.info("loaded {} samples:", Samples.size());
    }

    public boolean multipleSamples() { return Samples.size() > 1; }

    public String formFilename(final String fileType)
    {
        String fileName = OutputDir;

        if(multipleSamples())
        {
            fileName += "ctdna_cohort.";
        }
        else if(Samples.get(0).CtDnaSamples.size() > 1)
        {
            fileName += format("%s.ctdna.", Samples.get(0).PatientId);
        }
        else
        {
            fileName += format("%s_%s.ctdna.", Samples.get(0).PatientId, Samples.get(0).CtDnaSamples.get(0));
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

        StringJoiner sj = new StringJoiner(", ");
        Arrays.stream(PurityMethod.values()).forEach(x -> sj.add(x.toString()));
        configBuilder.addConfigItem(
                PURITY_METHODS, false,
                "List of purity methods separated by ',' default(all) from: " + sj);

        configBuilder.addConfigItem(SOMATIC_VCF, false, "Somatic VCF files, separated by ','", "");
        configBuilder.addConfigItem(SAMPLE_DATA_DIR_CFG, true, SAMPLE_DATA_DIR_DESC);
        configBuilder.addConfigItem(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        configBuilder.addConfigItem(COBALT_DIR_CFG, false, COBALT_DIR_DESC);
        configBuilder.addFlag(APPLY_DROPOUT, "Apply somatic drop-out logic");
        configBuilder.addFlag(WRITE_VARIANTS, "Write variants");
        configBuilder.addFlag(WRITE_CN_RATIOS, "Write copy number segment GC ratio summary");
        configBuilder.addFlag(PLOT_CN,"Plot copy number / GC ratio fit");
        configBuilder.addFlag(INCLUDE_FILTERED_VARIANTS, "Include filtered somatic variants in output (not purity calcs)");

        addRefGenomeConfig(configBuilder, false);

        configBuilder.addDecimal(
                NOISE_READS_PER_MILLION, "Expected reads-per-million from noise", DEFAULT_NOISE_READS_PER_MILLION);

        configBuilder.addDecimal(
                NOISE_READS_PER_MILLION_DUAL,
                "Expected reads-per-million from noise for dual-strand reads", DEFAULT_NOISE_READS_PER_MILLION_DUAL_STRAND);

        configBuilder.addDecimal( GC_RATIO_MIN,"GC ratio minimum permitted", DEFAULT_GC_RATIO_MIN);

        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
    }
}
