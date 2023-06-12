package com.hartwig.hmftools.ctdna.purity;

import static com.hartwig.hmftools.common.utils.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.DEFAULT_NOISE_READS_PER_MILLION;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.DEFAULT_NOISE_READS_PER_MILLION_DUAL_STRAND;
import static com.hartwig.hmftools.ctdna.purity.SampleData.ctDnaSamplesFromStr;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

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
    public final boolean WriteSomatics;
    public final boolean WriteCnRatios;
    public final boolean PlotCnFit;
    public final boolean WriteFilteredSomatics;
    public final double NoiseReadsPerMillion;
    public final double NoiseReadsPerMillionDualStrand;
    public final int Threads;

    private static final String PATIENT_ID = "patient_id";
    private static final String TUMOR_ID = "tumor_id";
    private static final String CTDNA_SAMPLES = "ctdna_samples";
    private static final String PURITY_METHODS = "purity_methods";
    private static final String SOMATIC_VCF = "somatic_vcf";
    private static final String SAMPLE_DATA_DIR = "sample_data_dir";
    private static final String PURPLE_DIR = "purple_dir";
    private static final String COBALT_DIR = "cobalt_dir";
    private static final String WRITE_VARIANTS = "write_somatics";
    private static final String INCLUDE_FILTERED_VARIANTS = "write_filtered_somatics";
    private static final String WRITE_CN_RATIOS = "write_cn_ratios";
    private static final String PLOT_CN = "plot_cn_fit";
    private static final String NOISE_READS_PER_MILLION = "noise_per_mill";
    private static final String NOISE_READS_PER_MILLION_DUAL = "noise_per_mill_dual";

    public PurityConfig(final CommandLine cmd)
    {
        SampleDataDir = checkAddDirSeparator(cmd.getOptionValue(SAMPLE_DATA_DIR));

        Samples = Lists.newArrayList();
        loadSampleData(cmd);

        PurityMethods = Lists.newArrayList();

        if(cmd.hasOption(PURITY_METHODS))
        {
            Arrays.stream(cmd.getOptionValue(PURITY_METHODS).split(ITEM_DELIM, -1))
                    .forEach(x -> PurityMethods.add(PurityMethod.valueOf(x)));
        }
        else
        {
            Arrays.stream(PurityMethod.values()).forEach(x -> PurityMethods.add(x));
        }

        SomaticVcf = cmd.getOptionValue(SOMATIC_VCF, "");
        PurpleDir = checkAddDirSeparator(cmd.getOptionValue(PURPLE_DIR, SampleDataDir));
        CobaltDir = checkAddDirSeparator(cmd.getOptionValue(COBALT_DIR, SampleDataDir));
        OutputDir = checkAddDirSeparator(cmd.getOptionValue(OUTPUT_DIR, SampleDataDir));
        OutputId = cmd.getOptionValue(OUTPUT_ID);

        NoiseReadsPerMillion = Double.parseDouble(
                cmd.getOptionValue(NOISE_READS_PER_MILLION, String.valueOf(DEFAULT_NOISE_READS_PER_MILLION)));
        NoiseReadsPerMillionDualStrand = Double.parseDouble(
                cmd.getOptionValue(NOISE_READS_PER_MILLION_DUAL, String.valueOf(DEFAULT_NOISE_READS_PER_MILLION_DUAL_STRAND)));

        WriteSomatics = cmd.hasOption(WRITE_VARIANTS);
        WriteCnRatios = cmd.hasOption(WRITE_CN_RATIOS);
        WriteFilteredSomatics = cmd.hasOption(INCLUDE_FILTERED_VARIANTS);
        PlotCnFit = cmd.hasOption(PLOT_CN);
        Threads = parseThreads(cmd);
    }

    private void loadSampleData(final CommandLine cmd)
    {
        if(cmd.hasOption(SAMPLE_ID_FILE))
        {
            String filename = cmd.getOptionValue(SAMPLE_ID_FILE);

            if(!Files.exists(Paths.get(filename)))
                filename = SampleDataDir + filename;

            try
            {
                final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

                String header = fileContents.get(0);
                Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, CSV_DELIM);
                fileContents.remove(0);

                int patientIndex = fieldsIndexMap.get("PatientId");
                int tumorIndex = fieldsIndexMap.get("TumorId");
                int ctdnaIndex = fieldsIndexMap.get("CtDnaSampleIds");
                Integer vcfIndex = fieldsIndexMap.get("VcfTag");

                for(String line : fileContents)
                {
                    if(line.startsWith("#") || line.isEmpty())
                        continue;

                    String[] values = line.split(CSV_DELIM, -1);

                    String vcfTag = vcfIndex != null && vcfIndex < values.length ? values[vcfIndex] : "";

                    Samples.add(new SampleData(
                            values[patientIndex], values[tumorIndex], ctDnaSamplesFromStr(values[ctdnaIndex]), vcfTag));
                }
            }
            catch (IOException e)
            {
                CT_LOGGER.error("failed to read sample data file({}): {}", filename, e.toString());
            }
        }
        else
        {
            Samples.add(new SampleData(
                    cmd.getOptionValue(PATIENT_ID),
                    cmd.getOptionValue(TUMOR_ID),
                    ctDnaSamplesFromStr(cmd.getOptionValue(CTDNA_SAMPLES)), ""));
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
        else
        {
            fileName += Samples.get(0).PatientId + ".ctdna.";
        }

        fileName += fileType;

        if(OutputId != null)
            fileName += "." + OutputId;

        fileName += TSV_EXTENSION;

        return fileName;
    }

    public static void addCommandLineOptions(final Options options)
    {
        options.addOption(SAMPLE_ID_FILE, true, "Patient and sample data file: PatientId,TumorId,CtDnaSampleIds");
        options.addOption(PATIENT_ID, true, "Patient ID");
        options.addOption(TUMOR_ID, true, "Original tumor sample ID");
        options.addOption(CTDNA_SAMPLES, true, "List of ctDNA sample IDs separated by ','");

        options.addOption(
                PURITY_METHODS, true,
                "List of purity methods separated by ',' default(all) from: " + PurityMethod.values());

        options.addOption(SOMATIC_VCF, true, "Somatic VCF files, separated by ','");
        options.addOption(SAMPLE_DATA_DIR, true, "Sample data directory for all files");
        options.addOption(PURPLE_DIR, true, "Sample Purple directory");
        options.addOption(COBALT_DIR, true, "Sample Cobalt directory");
        options.addOption(WRITE_VARIANTS, false, "Write variants");
        options.addOption(WRITE_CN_RATIOS, false, "Write copy number segment GC ratio summary");
        options.addOption(PLOT_CN, false, "Plot copy number / GC ratio fit");
        options.addOption(INCLUDE_FILTERED_VARIANTS, false, "Include filtered somatic variants in output (not purity calcs)");

        options.addOption(
                NOISE_READS_PER_MILLION, true,
                "Expected reads-per-million from noise, default: " + DEFAULT_NOISE_READS_PER_MILLION);

        addOutputOptions(options);
        addLoggingOptions(options);
        addThreadOptions(options);
    }
}
