package com.hartwig.hmftools.purple.tools;

import java.io.File;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;

public class PositiveAllelicDepth implements AutoCloseable
{
    private static final Logger LOGGER = LogManager.getLogger(PositiveAllelicDepth.class);

    private static final String IN_VCF = "in";
    private static final String OUT_VCF = "out";

    public static void main(String[] args) throws ParseException
    {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String inputFilePath = cmd.getOptionValue(IN_VCF);
        final String outputFilePath = cmd.getOptionValue(OUT_VCF);

        if(outputFilePath == null || inputFilePath == null)
        {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("PositiveAllelicDepth", options);
            System.exit(1);
        }

        try (PositiveAllelicDepth app = new PositiveAllelicDepth(inputFilePath, outputFilePath))
        {
            app.run();
        }
    }

    private final VCFFileReader fileReader;
    private final VariantContextWriter fileWriter;

    private PositiveAllelicDepth(final String inputVCF, final String outputVCF)
    {
        this.fileReader = new VCFFileReader(new File(inputVCF), false);
        this.fileWriter = new VariantContextWriterBuilder().setOutputFile(outputVCF).build();
    }

    private void run()
    {
        fileWriter.writeHeader(fileReader.getFileHeader());
        for(final VariantContext context : fileReader)
        {

            final VariantContextBuilder builder = new VariantContextBuilder(context);
            final List<Genotype> genotypes = Lists.newArrayList();

            for(int genotypeIndex = 0; genotypeIndex < context.getGenotypes().size(); genotypeIndex++)
            {
                Genotype genotype = context.getGenotype(genotypeIndex);

                if(genotype.hasAD() && genotype.hasDP())
                {
                    int[] ad = genotype.getAD();
                    int total = 0;
                    boolean updateRecord = false;
                    for(int i = 0; i < ad.length; i++)
                    {
                        int adValue = ad[i];
                        if(adValue < 0)
                        {
                            updateRecord = true;
                            ad[i] = Math.abs(adValue);
                        }
                        total += Math.abs(adValue);
                    }

                    if(updateRecord)
                    {
                        LOGGER.info("Updated entry at {}:{}", context.getContig(), context.getStart());
                        Genotype updated = new GenotypeBuilder(genotype).AD(ad).DP(total).make();
                        genotypes.add(updated);
                    }
                    else
                    {
                        genotypes.add(genotype);
                    }
                }
            }

            builder.genotypes(genotypes);
            fileWriter.add(builder.make());
        }
    }

    @Override
    public void close()
    {
        fileReader.close();
        fileWriter.close();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(IN_VCF, true, "Input file.");
        options.addOption(OUT_VCF, true, "Output file.");
        return options;
    }
}
