package com.hartwig.hmftools.peach;

import static com.hartwig.hmftools.peach.TestUtils.getTestResourcePath;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PeachConfigTest
{
    @Test
    public void testMinimalCommandLineArguments()
    {
        String vcfFile = getTestResourcePath("variants.vcf.gz");
        String haplotypesFile = getTestResourcePath("haplotypes.complicated.37.tsv");
        String sampleName = "FAKENAME";
        String outputDir = "/path/to/output";

        String[] args = {
                "-vcf_file", vcfFile,
                "-haplotypes_file", haplotypesFile,
                "-sample_name", sampleName,
                "-output_dir", outputDir
        };

        PeachConfig config = constructPeachConfigFromArgs(args);

        assertEquals(vcfFile, config.vcfFile);
        assertEquals(haplotypesFile, config.haplotypesFile);
        assertEquals(sampleName, config.sampleName);
        assertEquals(outputDir + "/", config.outputDir);

        assertEquals("/path/to/output/FAKENAME.peach.events.tsv", config.getEventsOutputPath());
        assertEquals("/path/to/output/FAKENAME.peach.gene.events.tsv", config.getEventsPerGeneOutputPath());
        assertEquals("/path/to/output/FAKENAME.peach.haplotypes.all.tsv", config.getAllHaplotypeCombinationsOutputPath());
        assertEquals("/path/to/output/FAKENAME.peach.haplotypes.best.tsv", config.getBestHaplotypeCombinationsOutputPath());
        assertEquals("/path/to/output/FAKENAME.peach.qc.tsv", config.getQcStatusOutputPath());
    }

    @Test
    public void testMaximalCommandLineArguments()
    {
        String vcfFile = getTestResourcePath("variants.vcf.gz");
        String haplotypesFile = getTestResourcePath("haplotypes.complicated.37.tsv");
        String drugsFile = getTestResourcePath("drugs.tsv");
        String functionFile = getTestResourcePath("function.tsv");
        String sampleName = "FAKENAME";
        String outputDir = "/path/to/output";

        String[] args = {
                "-vcf_file", vcfFile,
                "-haplotypes_file", haplotypesFile,
                "-sample_name", sampleName,
                "-output_dir", outputDir,
                "-drugs_file", drugsFile,
                "-function_file", functionFile
        };

        PeachConfig config = constructPeachConfigFromArgs(args);

        assertEquals(vcfFile, config.vcfFile);
        assertEquals(haplotypesFile, config.haplotypesFile);
        assertEquals(sampleName, config.sampleName);
        assertEquals(outputDir + "/", config.outputDir);
        assertEquals(drugsFile, config.drugsFile);
        assertEquals(functionFile, config.functionFile);

        assertEquals("/path/to/output/FAKENAME.peach.events.tsv", config.getEventsOutputPath());
        assertEquals("/path/to/output/FAKENAME.peach.gene.events.tsv", config.getEventsPerGeneOutputPath());
        assertEquals("/path/to/output/FAKENAME.peach.haplotypes.all.tsv", config.getAllHaplotypeCombinationsOutputPath());
        assertEquals("/path/to/output/FAKENAME.peach.haplotypes.best.tsv", config.getBestHaplotypeCombinationsOutputPath());
        assertEquals("/path/to/output/FAKENAME.peach.qc.tsv", config.getQcStatusOutputPath());
    }

    @NotNull
    private static PeachConfig constructPeachConfigFromArgs(@NotNull String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder("Peach");
        PeachConfig.addOptions(configBuilder);
        configBuilder.checkAndParseCommandLine(args);
        return new PeachConfig(configBuilder);
    }
}
