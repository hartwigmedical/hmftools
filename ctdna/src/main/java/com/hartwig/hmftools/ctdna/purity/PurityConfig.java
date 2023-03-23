package com.hartwig.hmftools.ctdna.purity;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.DELIMETER;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.DEFAULT_NOISE_READS_PER_MILLION;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class PurityConfig
{
    public final String PatientId;
    public final String TumorId;
    public final List<String> CtDnaSamples;
    public final List<String> SomaticVcfs;
    public final String PurpleDir;
    public final String CobaltDir;
    public final String OutputDir;
    public final String OutputId;
    public final boolean WriteVariants;
    public final boolean WriteCnRatios;
    public final int NoiseReadsPerMillion;

    private static final String PATIENT_ID = "patient_id";
    private static final String TUMOR_ID = "tumor_id";
    private static final String CTDNA_SAMPLES = "ctdna_samples";
    private static final String SOMATIC_VCFS = "somatic_vcfs";
    private static final String PURPLE_DIR = "purple_dir";
    private static final String COBALT_DIR = "cobalt_dir";
    private static final String WRITE_VARIANTS = "write_variants";
    private static final String WRITE_CN_RATIOS = "write_cn_ratios";
    private static final String NOISE_READS_PER_MILLION = "noise_per_mill";

    public PurityConfig(final CommandLine cmd)
    {
        PatientId = cmd.getOptionValue(PATIENT_ID);
        TumorId = cmd.getOptionValue(TUMOR_ID);

        CT_LOGGER.info("patient({}) tumor({})", PatientId, TumorId);

        CtDnaSamples = Lists.newArrayList();

        String ctDnaSamples = cmd.getOptionValue(CTDNA_SAMPLES);
        String sampleDelim = ctDnaSamples.contains(DELIMETER) ? DELIMETER : ";";

        for(String ctDnaSample : ctDnaSamples.split(sampleDelim, -1))
        {
            CT_LOGGER.info("added ctDNA sample({})", ctDnaSample);
            CtDnaSamples.add(ctDnaSample);
        }

        SomaticVcfs = Lists.newArrayList();
        for(String vcf : cmd.getOptionValue(SOMATIC_VCFS).split(",", -1))
        {
            CT_LOGGER.info("added somatic VCF({})", vcf);
            SomaticVcfs.add(vcf);
        }

        PurpleDir = checkAddDirSeparator(cmd.getOptionValue(PURPLE_DIR));
        CobaltDir = checkAddDirSeparator(cmd.getOptionValue(COBALT_DIR));
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);

        NoiseReadsPerMillion = Integer.parseInt(cmd.getOptionValue(NOISE_READS_PER_MILLION, String.valueOf(DEFAULT_NOISE_READS_PER_MILLION)));

        WriteVariants = cmd.hasOption(WRITE_VARIANTS);
        WriteCnRatios = cmd.hasOption(WRITE_CN_RATIOS);
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(PATIENT_ID, true, "Patient ID");
        options.addOption(TUMOR_ID, true, "Original tumor sample ID");
        options.addOption(CTDNA_SAMPLES, true, "List of ctDNA sample IDs separated by ','");
        options.addOption(SOMATIC_VCFS, true, "Somatic VCF files, separated by ','");
        options.addOption(PURPLE_DIR, true, "Sample Purple directory");
        options.addOption(COBALT_DIR, true, "Sample Cobalt directory");
        options.addOption(WRITE_VARIANTS, false, "Write variants");
        options.addOption(WRITE_CN_RATIOS, false, "Write copy number segment GC ratio summary");

        options.addOption(
                NOISE_READS_PER_MILLION, true,
                "Expected reads-per-million from noise, default: " + DEFAULT_NOISE_READS_PER_MILLION);

        addOutputOptions(options);
        addLoggingOptions(options);
    }
}
