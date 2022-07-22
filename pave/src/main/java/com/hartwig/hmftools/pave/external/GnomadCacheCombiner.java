package com.hartwig.hmftools.pave.external;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class GnomadCacheCombiner
{
    private final String mGnomadDir1;
    private final String mOutputDir;
    private final String mGnomadDir2;

    private static final String GNOMAD_DIR_1 = "gnomad_dir_1";
    private static final String GNOMAD_DIR_2 = "gnomad_dir_2";

    public GnomadCacheCombiner(final CommandLine cmd)
    {
        mOutputDir = parseOutputDir(cmd);
        mGnomadDir1 = checkAddDirSeparator(cmd.getOptionValue(GNOMAD_DIR_1));
        mGnomadDir2 = checkAddDirSeparator(cmd.getOptionValue(GNOMAD_DIR_2));
    }

    public void run()
    {
        if(mGnomadDir1 == null || !Files.exists(Paths.get(mGnomadDir1)))
        {
            PV_LOGGER.error("missing input file, exiting");
            System.exit(1);
        }

        PV_LOGGER.info("combining Gnomad files from {} and {}, writing to {}", mGnomadDir1, mGnomadDir2, mOutputDir);

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            combineChromosome(chromosome.toString());
        }

        PV_LOGGER.info("Gnomad file combine complete");
    }

    private void combineChromosome(final String chromosome)
    {
        // eg gnomad_variants_chr21_v38.csv.gz
        String chrFilename = format("gnomad_variants_chr%s_v38.csv.gz", chromosome);

        String filename1 = mGnomadDir1 + chrFilename;
        String filename2 = mGnomadDir2 + chrFilename;

        if(!Files.exists(Paths.get(filename1)) || !Files.exists(Paths.get(filename2)))
        {
            PV_LOGGER.info("skipping chromosome({})", chromosome);
            return;
        }

        PV_LOGGER.info("combining chromosome({})", chromosome);

        try
        {

            BufferedReader fileReader = createBufferedReader(filename1);

            int itemCount = 0;
            String line = fileReader.readLine(); // skip header

            List<GnomadVariant> variantsList = Lists.newArrayList();

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(",", -1);

                int position = Integer.parseInt(values[0]);
                String ref = values[1];
                String alt = values[2];
                double frequency = Double.parseDouble(values[3]);

                variantsList.add(new GnomadVariant(position, ref, alt, frequency));
                ++itemCount;
            }

            PV_LOGGER.info("loaded {} Gnomad frequency records from file({})", itemCount, filename1);

            // load and combine the second

            fileReader = createBufferedReader(filename2);

            itemCount = 0;
            int matchedCount = 0;

            line = fileReader.readLine(); // skip header
            int currentIndex = 0;

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(",", -1);

                ++itemCount;

                int position = Integer.parseInt(values[0]);
                String ref = values[1];
                String alt = values[2];
                double frequency = Double.parseDouble(values[3]);

                boolean matched = false;

                while(currentIndex < variantsList.size())
                {
                    GnomadVariant variant = variantsList.get(currentIndex);

                    if(variant.Position < position)
                    {
                        ++currentIndex;
                        continue;
                    }
                    else if(variant.Position == position)
                    {
                        if(variant.matches(position, ref, alt))
                        {
                            matched = true;
                            ++matchedCount;
                            variant.Frequency = max(variant.Frequency, frequency);
                            break;
                        }

                        ++currentIndex;
                        continue;
                    }

                    if(variantsList.get(currentIndex).Position > position)
                        break;
                }

                if(!matched)
                    variantsList.add(currentIndex, new GnomadVariant(position, ref, alt, frequency));

                if((itemCount % 100000) == 0)
                {
                    PV_LOGGER.debug("chr({}) {} items combined", chromosome, itemCount);
                }
            }

            PV_LOGGER.info("loaded {} Gnomad frequency records from file({}), matched({})", itemCount, filename1, matchedCount);

            String outputFile = mOutputDir + chrFilename;
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("Position,Ref,Alt,Frequency");
            writer.newLine();

            for(GnomadVariant variant : variantsList)
            {
                writer.write(String.format("%d,%s,%s,%.5f", variant.Position, variant.Ref, variant.Alt, variant.Frequency));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to load Gnomad frequency file: {}", e.toString());
            return;
        }
    }

    private class GnomadVariant
    {
        public final int Position;
        public final String Ref;
        public final String Alt;
        public double Frequency;

        public GnomadVariant(final int position, final String ref, final String alt, final double frequency)
        {
            Position = position;
            Ref = ref;
            Alt = alt;
            Frequency = frequency;
        }

        public boolean matches(final int position, final String ref, final String alt)
        {
            return Position == position && ref.equals(Ref) && alt.equals(Alt);
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        options.addOption(GNOMAD_DIR_1, true, "Gnomad VCF input file");
        options.addOption(GNOMAD_DIR_2, true, "Gnomad VCF input file");
        addOutputOptions(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);
        setLogLevel(cmd);

        GnomadCacheCombiner gnomadCacheBuilder = new GnomadCacheCombiner(cmd);
        gnomadCacheBuilder.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
