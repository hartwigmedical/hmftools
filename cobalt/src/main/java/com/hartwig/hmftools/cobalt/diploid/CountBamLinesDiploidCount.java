package com.hartwig.hmftools.cobalt.diploid;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.MAX_DIPLOID;
import static com.hartwig.hmftools.cobalt.CobaltConstants.MIN_DIPLOID;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatioFactory;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class CountBamLinesDiploidCount implements AutoCloseable
{
    private final String mInputFile;
    private final String mOutputFile;
    private final long mTimestamp = System.currentTimeMillis();

    public CountBamLinesDiploidCount(final String inputFile, final String outputFile)
    {
        mInputFile = inputFile;
        mOutputFile = outputFile;
    }

    public void run() throws IOException
    {
        boolean outputFileExists = new File(mOutputFile).exists();
        final Map<GenomePosition, DiploidCount> diploidCountMap;
        if(outputFileExists)
        {
            CB_LOGGER.info("Reading existing output: {}", mOutputFile);
            diploidCountMap = DiploidCount.readDiploidCountAsMap(mOutputFile);
        }
        else
        {
            diploidCountMap = Maps.newHashMap();
        }

        CB_LOGGER.info("Reading input file {}", mInputFile);
        final ListMultimap<Chromosome, CobaltRatio> ratios = CobaltRatioFile.read(mInputFile);
        final List<MedianRatio> medianRatios = MedianRatioFactory.create(ratios);
        final CobaltChromosomes chromosomes = new CobaltChromosomes(medianRatios);
        if(chromosomes.hasGermlineAberrations())
        {
            CB_LOGGER.warn("Excluding {} because of germline aberrations", mInputFile);
            return;
        }

        Collection<CobaltChromosome> selectedChromosomes;
        if(chromosomes.contains("Y"))
        {
            selectedChromosomes = Collections.singleton(chromosomes.get("Y"));
        }
        else if(chromosomes.contains("chrY"))
        {
            selectedChromosomes = Collections.singleton(chromosomes.get("chrY"));
        }
        else
        {
            selectedChromosomes = chromosomes.chromosomes();
        }

        CB_LOGGER.info("Processing diploid regions");
        for(CobaltChromosome selectedChromosome : selectedChromosomes)
        {
            Collection<CobaltRatio> selectedRatios = ratios.get(HumanChromosome.fromString(selectedChromosome.contig()));
            for(CobaltRatio ratio : selectedRatios)
            {
                GenomePosition position = GenomePositions.create(ratio.chromosome(), ratio.position());
                DiploidCount count = diploidCountMap.computeIfAbsent(position, x -> new DiploidCount(x, 0, 0));
                count.incrementTotal();
                if(isDiploid(selectedChromosome, ratio))
                {
                    count.incrementDiploid();
                }
            }
        }

        CB_LOGGER.info("Writing output: {}", mOutputFile);
        Files.write(new File(mOutputFile).toPath(),
                diploidCountMap.values().stream().sorted().map(DiploidCount::toString).collect(Collectors.toList()));
    }

    private static boolean isDiploid(final CobaltChromosome chromosome, CobaltRatio ratio)
    {
        double value = ratio.referenceGCDiploidRatio();
        return (Doubles.greaterOrEqual(value, MIN_DIPLOID * chromosome.actualRatio()) && Doubles.lessOrEqual(value,
                MAX_DIPLOID * chromosome.actualRatio()));
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

        try(CountBamLinesDiploidCount app = new CountBamLinesDiploidCount(inputFile, outputFile))
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
        options.addOption("in", true, "Input ratios");
        options.addOption("out", true, "Output count");
        return options;
    }
}
