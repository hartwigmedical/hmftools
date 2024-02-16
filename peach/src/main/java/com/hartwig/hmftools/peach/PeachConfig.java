package com.hartwig.hmftools.peach;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import org.jetbrains.annotations.NotNull;

public class PeachConfig
{
    public final String vcfFile;
    public final String sampleName;
    public final String haplotypesFile;
    public final String drugsFile;
    public final String functionFile;
    public final String outputDir;
    private static final String VCF_FILE = "vcf_file";
    private static final String SAMPLE_NAME = "sample_name";
    private static final String HAPLOTYPES_FILE = "haplotypes_file";
    private static final String DRUGS_FILE = "drugs_file";
    private static final String FUNCTION_FILE = "function_file";

    public PeachConfig(@NotNull ConfigBuilder configBuilder)
    {
        vcfFile = configBuilder.getValue(VCF_FILE);
        haplotypesFile = configBuilder.getValue(HAPLOTYPES_FILE);
        sampleName = configBuilder.getValue(SAMPLE_NAME);
        drugsFile = configBuilder.hasValue(DRUGS_FILE) ? configBuilder.getValue(DRUGS_FILE) : null;
        functionFile = configBuilder.hasValue(FUNCTION_FILE) ? configBuilder.getValue(FUNCTION_FILE) : null;

        outputDir = parseOutputDir(configBuilder);
    }

    public boolean isValid()
    {
        return vcfFile != null && sampleName != null && haplotypesFile != null && outputDir != null;
    }

    public static void addOptions(@NotNull ConfigBuilder configBuilder)
    {
        configBuilder.addPath(VCF_FILE, true, "VCF input file");
        configBuilder.addPath(HAPLOTYPES_FILE, true, "Haplotype config file");
        configBuilder.addConfigItem(SAMPLE_NAME, true, "Name of sample in VCF to call haplotypes for");
        configBuilder.addPath(DRUGS_FILE, false, "Config file of relevant drugs");
        configBuilder.addPath(FUNCTION_FILE, false, "Config file for haplotype function");
        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);
    }

    @NotNull
    public String getEventsOutputPath()
    {
        return outputDir + sampleName + ".peach.events.tsv";
    }

    @NotNull
    public String getEventsPerGeneOutputPath()
    {
        return outputDir + sampleName + ".peach.gene.events.tsv";
    }

    @NotNull
    public String getAllHaplotypeCombinationsOutputPath()
    {
        return outputDir + sampleName + ".peach.haplotypes.all.tsv";
    }

    @NotNull
    public String getBestHaplotypeCombinationsOutputPath()
    {
        return outputDir + sampleName + ".peach.haplotypes.best.tsv";
    }

    @NotNull
    public String getQcStatusOutputPath()
    {
        return outputDir + sampleName + ".peach.qc.tsv";
    }
}
