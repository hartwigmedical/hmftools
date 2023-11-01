package com.hartwig.hmftools.peach;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import java.io.File;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

public class PeachConfig
{
    public final String vcfFile;
    public final String sampleName;
    public final String haplotypesTsv;
    public final String outputDir;

    public final boolean doLiftOver;
    public final String chainFile;
    public final String liftOverBed;
    public final String picardJar;
    public final String targetRefGenome;

    private static final String VCF_FILE = "vcf_file";
    private static final String SAMPLE_NAME = "sample_name";
    private static final String HAPLOTYPES_TSV = "haplotypes_tsv";

    private static final String DO_LIFT_OVER = "do_lift_over";
    private static final String CHAIN_FILE = "chain_file";
    private static final String LIFT_OVER_BED = "lift_over_bed";
    private static final String PICARD = "picard_jar";
    private static final String TARGET_REF_GENOME = "target_ref_genome";

    public PeachConfig(final ConfigBuilder configBuilder)
    {
        vcfFile = configBuilder.getValue(VCF_FILE);
        haplotypesTsv = configBuilder.getValue(HAPLOTYPES_TSV);
        sampleName = configBuilder.getValue(SAMPLE_NAME);

        doLiftOver = configBuilder.hasFlag(DO_LIFT_OVER);
        chainFile = configBuilder.getValue(CHAIN_FILE);
        liftOverBed = configBuilder.getValue(LIFT_OVER_BED);
        picardJar = configBuilder.getValue(PICARD);
        targetRefGenome = configBuilder.getValue(TARGET_REF_GENOME);

        outputDir = parseOutputDir(configBuilder);
    }

    public boolean isValid()
    {
        if (vcfFile == null || sampleName == null || haplotypesTsv == null || outputDir == null)
            return false;
        if (doLiftOver)
            return !(liftOverBed == null || picardJar == null || chainFile == null || targetRefGenome == null);

        return true;
    }

    public static void addOptions(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(VCF_FILE, true, "VCF input file");
        configBuilder.addPath(HAPLOTYPES_TSV, true, "Haplotype config file");
        configBuilder.addConfigItem(SAMPLE_NAME, true, "Name of sample in VCF to call haplotypes for");
        configBuilder.addFlag(DO_LIFT_OVER, "Do liftover to 38");
        configBuilder.addPath(CHAIN_FILE, false, "USCS chain file for liftover, if liftover is needed");
        configBuilder.addPath(LIFT_OVER_BED, false, "BED file specifying ranges that need to be lifted over to 38, if liftover is needed");
        configBuilder.addPath(PICARD, false, "Picard JAR for liftover, if liftover is needed");
        configBuilder.addPath(TARGET_REF_GENOME, false, "Target reference FASTA file for liftover, if liftover is needed");

        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);
    }

    public String getAdjustedChainFilePath()
    {
        return getExtendedFileName(outputDir, chainFile, "adjusted", ".over");
    }

    public String getLiftoverOutputVcfPath()
    {
        return getExtendedFileName(outputDir, vcfFile, "liftover", ".vcf");
    }

    public String getLiftoverRejectVcfPath()
    {
        return getExtendedFileName(outputDir, vcfFile, "liftover_reject", ".vcf");
    }

    public String getEventsOutputPath()
    {
        return outputDir + sampleName + ".peach.events.tsv";
    }

    public String getEventsPerGeneOutputPath()
    {
        return outputDir + sampleName + ".peach.gene.events.tsv";
    }

    public String getAllHaplotypeCombinationsOutputPath()
    {
        return outputDir + sampleName + ".peach.haplotypes.all.tsv";
    }

    public String getBestHaplotypeCombinationsOutputPath()
    {
        return outputDir + sampleName + ".peach.haplotypes.best.tsv";
    }

    public String getQcStatusOutputPath()
    {
        return outputDir + sampleName + ".peach.qc.tsv";
    }

    public static String getExtendedFileName(String outputDir, String originalFileName, String addition, String addBefore)
    {
        File file = new File(originalFileName);
        String filename = file.getName();
        int extensionIndex = filename.indexOf(addBefore);
        return outputDir + filename.substring(0, extensionIndex) + "." + addition + filename.substring(extensionIndex);
    }
}
