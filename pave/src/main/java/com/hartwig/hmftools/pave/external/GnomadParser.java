package com.hartwig.hmftools.pave.external;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Collectors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class GnomadParser
{
    private final String mInputVcf;
    private final String mOutputDir;

    private static final String GNOMAD_FILE = "gnomad_file";

    public GnomadParser(final CommandLine cmd)
    {
        mOutputDir = parseOutputDir(cmd);
        mInputVcf = cmd.getOptionValue(GNOMAD_FILE);
    }

    public void run()
    {
        if(mInputVcf == null || !Files.exists(Paths.get(mInputVcf)))
        {
            PV_LOGGER.error("missing input file, exiting");
            System.exit(1);
        }

        PV_LOGGER.info("parsing Gnomad file({})", mInputVcf);

        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(
                    mInputVcf, new VCFCodec(), false);

            String outputFile = mOutputDir + "gnomad_variants.csv";
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("Chromosome,Position,Ref,Alt,Frequency");
            writer.newLine();

            int itemCount = 0;
            int errorCount = 0;

            for(VariantContext context : reader.iterator())
            {
                try
                {
                    String chromosome = context.getContig();
                    int position = context.getStart();

                    String ref = context.getReference().getBaseString();
                    String alt = context.getAlternateAlleles().stream().map(Allele::toString).collect(Collectors.joining(","));

                    double frequency = context.getAttributeAsDouble("AF", 0);

                    writer.write(String.format("%s,%d,%s,%s,%4.5e", chromosome, position, ref, alt, frequency));
                    writer.newLine();
                }
                catch(Exception e)
                {
                    ++errorCount;
                    PV_LOGGER.debug("error processing variant at item({})", itemCount);

                    if(errorCount > 100)
                        break;
                }

                ++itemCount;

                if(itemCount > 0 && (itemCount % 10000) == 0)
                {
                    PV_LOGGER.debug("processed {} variants", itemCount);
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            PV_LOGGER.error(" failed to read gnomad file): {}", e.toString());
        }

        PV_LOGGER.info("Gnomad file parse complete");
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        options.addOption(GNOMAD_FILE, true, "Gnomad VCF input file");
        addOutputDir(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);
        setLogLevel(cmd);

        GnomadParser gnomadParser = new GnomadParser(cmd);
        gnomadParser.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
