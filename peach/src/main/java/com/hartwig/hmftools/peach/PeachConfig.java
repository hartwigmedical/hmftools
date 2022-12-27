package com.hartwig.hmftools.peach;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

import java.io.File;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

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

    public PeachConfig(final CommandLine cmd)
    {
        vcfFile = cmd.getOptionValue(VCF_FILE);
        haplotypesTsv = cmd.getOptionValue(HAPLOTYPES_TSV);
        sampleName = cmd.getOptionValue(SAMPLE_NAME);

        doLiftOver = cmd.hasOption(DO_LIFT_OVER);
        chainFile = cmd.getOptionValue(CHAIN_FILE);
        liftOverBed = cmd.getOptionValue(LIFT_OVER_BED);
        picardJar = cmd.getOptionValue(PICARD);
        targetRefGenome = cmd.getOptionValue(TARGET_REF_GENOME);

        outputDir = parseOutputDir(cmd);
    }

    public boolean isValid()
    {
        if (vcfFile == null || sampleName == null || haplotypesTsv == null || outputDir == null)
            return false;
        if (doLiftOver)
            return !(liftOverBed == null || picardJar == null || chainFile == null || targetRefGenome == null);

        return true;
    }

    public static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(VCF_FILE, true, "VCF input file");
        options.addOption(HAPLOTYPES_TSV, true, "Haplotype config file");
        options.addOption(SAMPLE_NAME, true, "Name of sample in VCF to call haplotypes for");

        options.addOption(DO_LIFT_OVER, false, "Do liftover to 38");
        options.addOption(CHAIN_FILE, true, "USCS chain file for liftover, if liftover is needed");
        options.addOption(LIFT_OVER_BED, true, "BED file specifying ranges that need to be lifted over to 38, if liftover is needed");
        options.addOption(PICARD, true, "Picard JAR for liftover, if liftover is needed");
        options.addOption(TARGET_REF_GENOME, true, "Target reference FASTA file for liftover, if liftover is needed");

        addOutputDir(options);
        addLoggingOptions(options);

        return options;
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

    public static String getExtendedFileName(String outputDir, String originalFileName, String addition, String addBefore)
    {
        File file = new File(originalFileName);
        String filename = file.getName();
        int extensionIndex = filename.indexOf(addBefore);
        return outputDir + filename.substring(0, extensionIndex) + "." + addition + filename.substring(extensionIndex);
    }
}
