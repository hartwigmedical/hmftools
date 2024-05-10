package com.hartwig.hmftools.peach;

import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PeachConfig
{
    @NotNull
    public final String vcfFile;
    @NotNull
    public final String sampleName;
    @NotNull
    public final String haplotypesFile;
    @NotNull
    public final String outputDir;
    @Nullable
    public final String drugsFile;
    @Nullable
    public final String functionFile;

    private static final String VCF_FILE = "vcf_file";
    private static final String SAMPLE_NAME = "sample_name";
    private static final String HAPLOTYPES_FILE = "haplotypes_file";
    private static final String DRUGS_FILE = "drugs_file";
    private static final String FUNCTION_FILE = "function_file";

    public PeachConfig(@NotNull ConfigBuilder configBuilder)
    {
        String nullableVcfFile = configBuilder.getValue(VCF_FILE);
        String nullableHaplotypesFile = configBuilder.getValue(HAPLOTYPES_FILE);
        String nullableSampleName = configBuilder.getValue(SAMPLE_NAME);
        String nullableOutputDir = parseOutputDir(configBuilder);

        if(nullableVcfFile == null || nullableHaplotypesFile == null || nullableSampleName == null || nullableOutputDir == null)
        {
            throw new IllegalArgumentException("invalid config");
        }

        vcfFile = nullableVcfFile;
        haplotypesFile = nullableHaplotypesFile;
        sampleName = nullableSampleName;
        outputDir = nullableOutputDir;

        drugsFile = configBuilder.hasValue(DRUGS_FILE) ? configBuilder.getValue(DRUGS_FILE) : null;
        functionFile = configBuilder.hasValue(FUNCTION_FILE) ? configBuilder.getValue(FUNCTION_FILE) : null;
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
        return checkAddDirSeparator(outputDir) + sampleName + ".peach.events.tsv";
    }

    @NotNull
    public String getEventsPerGeneOutputPath()
    {
        return checkAddDirSeparator(outputDir) + sampleName + ".peach.gene.events.tsv";
    }

    @NotNull
    public String getAllHaplotypeCombinationsOutputPath()
    {
        return checkAddDirSeparator(outputDir) + sampleName + ".peach.haplotypes.all.tsv";
    }

    @NotNull
    public String getBestHaplotypeCombinationsOutputPath()
    {
        return PeachGenotypeFile.generateFileName(outputDir, sampleName);
    }

    @NotNull
    public String getQcStatusOutputPath()
    {
        return checkAddDirSeparator(outputDir) + sampleName + ".peach.qc.tsv";
    }
}
