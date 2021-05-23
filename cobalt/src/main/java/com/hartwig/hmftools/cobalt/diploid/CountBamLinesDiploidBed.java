package com.hartwig.hmftools.cobalt.diploid;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class CountBamLinesDiploidBed implements AutoCloseable
{
    private final String mInputFile;
    private final String mOutputFile;
    private final long mTimestamp = System.currentTimeMillis();

    public CountBamLinesDiploidBed(final String inputFile, final String outputFile)
    {
        mInputFile = inputFile;
        mOutputFile = outputFile;
    }

    public void run() throws IOException
    {
        double cutoff = 0.50;

        CB_LOGGER.info("Reading input file: {}", mInputFile);
        final List<DiploidCount> diploidCounts = DiploidCount.readDiploidCountAsList(mInputFile);

        CB_LOGGER.info("Determining diploid regions with {} cutoff", cutoff);

        final DiploidRegionBuilder builder = new DiploidRegionBuilder(cutoff, 62, 34);
        diploidCounts.forEach(builder);

        final List<GenomeRegion> diploidRegions = builder.build();
        CB_LOGGER.info("Total diploid bases: {}", builder.getTotalDiploidBases());

        Files.write(new File(mOutputFile).toPath(),
                diploidRegions.stream().sorted().map(CountBamLinesDiploidBed::toBedString).collect(Collectors.toList()));
    }

    private static String toBedString(final GenomeRegion region)
    {
        return new StringJoiner(DELIMITER)
                .add(region.chromosome())
                .add(String.valueOf(region.start() - 1))
                .add(String.valueOf(region.end()))
                .toString();
    }

    @Override
    public void close()
    {
        CB_LOGGER.info("Complete in {} seconds", (System.currentTimeMillis() - mTimestamp) / 1000);
    }

    public static void main(String[] args) throws ParseException
    {
        final Options options = createOptions();
        final CommandLine cmd = new DefaultParser().parse(options, args);
        final String inputFile = cmd.getOptionValue("in");
        final String outputFile = cmd.getOptionValue("out");

        try(CountBamLinesDiploidBed app = new CountBamLinesDiploidBed(inputFile, outputFile))
        {
            app.run();
        }
        catch(Exception e)
        {
            CB_LOGGER.error(e);
            e.printStackTrace();
            System.exit(-1);
        }
    }

    @NotNull
    private static Options createOptions()
    {
        final Options options = new Options();
        options.addOption("in", true, "Cobalt diploid count");
        options.addOption("out", true, "Bed File");
        return options;
    }
}
