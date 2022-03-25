package com.hartwig.hmftools.common.utils.pcf;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Comparator;
import java.util.List;
import java.util.function.Function;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.genome.window.Window;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PCFFile
{
    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "sampleID";
    private static final String RATIO_EXTENSION = ".cobalt.ratio.pcf";
    private static final String BAF_EXTENSION = ".amber.baf.pcf";

    @NotNull
    public static String generateRatioFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + RATIO_EXTENSION;
    }

    @NotNull
    public static String generateBAFFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + BAF_EXTENSION;
    }

    @NotNull
    public static ListMultimap<Chromosome, PCFPosition> readPositions(int windowSize, final PCFSource source, final String filename) throws IOException
    {
        ListMultimap<Chromosome, PCFPosition> result = ArrayListMultimap.create();
        final Window window = new Window(windowSize);

        PCFPosition pcfPosition = null;

        String prevChromosome = Strings.EMPTY;
        int minPosition = 1;
        List<PCFPosition> chromosomeResult = Lists.newArrayList();

        for(String line : Files.readAllLines(new File(filename).toPath()))
        {
            if(!line.startsWith(HEADER_PREFIX))
            {
                String[] values = line.split(DELIMITER);
                final String chromosomeName = values[1];
                if(HumanChromosome.contains(chromosomeName))
                {
                    if(!chromosomeName.equals(prevChromosome))
                    {
                        if(pcfPosition != null)
                        {
                            chromosomeResult.add(pcfPosition);
                            result.putAll(HumanChromosome.fromString(prevChromosome), mergePositions(chromosomeResult));
                        }
                        chromosomeResult.clear();
                        pcfPosition = null;
                        minPosition = 1;
                        prevChromosome = chromosomeName;
                    }

                    int start = window.start(Integer.parseInt(values[3]));
                    int end = window.start(Integer.parseInt(values[4])) + windowSize;
                    if(pcfPosition != null)
                    {
                        pcfPosition.setMaxPosition(start);
                        chromosomeResult.add(pcfPosition);
                    }

                    pcfPosition = new PCFPosition(source, chromosomeName, start);
                    pcfPosition.setMinPosition(minPosition);
                    pcfPosition.setMaxPosition(start);

                    chromosomeResult.add(pcfPosition);

                    minPosition = end;

                    pcfPosition = new PCFPosition(source, chromosomeName, end);
                    pcfPosition.setMinPosition(end);
                    pcfPosition.setMaxPosition(end);
                }
            }
        }

        if(pcfPosition != null)
        {
            chromosomeResult.add(pcfPosition);
            result.putAll(HumanChromosome.fromString(prevChromosome), mergePositions(chromosomeResult));
        }

        return result;
    }

    public static Multimap<String, GenomeRegion> read(int windowSize, final String filename) throws IOException
    {
        return fromLines(windowSize, Files.readAllLines(new File(filename).toPath()));
    }

    private static Multimap<String, GenomeRegion> fromLines(int windowSize, List<String> lines)
    {
        Multimap<String, GenomeRegion> result = ArrayListMultimap.create();
        for(String line : lines)
        {
            if(!line.startsWith(HEADER_PREFIX))
            {
                final GenomeRegion region = fromString(windowSize, line);
                result.put(region.chromosome(), region);
            }
        }
        return result;
    }

    private static GenomeRegion fromString(int windowSize, final String line)
    {
        String[] values = line.split(DELIMITER);
        String chromosome = values[1];
        int posStart = Integer.parseInt(values[3]);
        int posEnd = Integer.parseInt(values[4]) + windowSize - 1;
        return GenomeRegions.create(chromosome, posStart, posEnd);
    }

    private static final Comparator<PCFPosition> COMPARE = Comparator.comparing(
            (Function<PCFPosition, String>) GenomePosition::chromosome).thenComparingLong(GenomePosition::position);

    private static List<PCFPosition> mergePositions(final List<PCFPosition> positions)
    {
        positions.sort(COMPARE);

        int i = 0;
        while(i < positions.size() - 1)
        {
            PCFPosition current = positions.get(i);
            PCFPosition next = positions.get(i + 1);

            if(current.chromosome().equals(next.chromosome()) && current.position() == next.position())
            {
                current.setMinPosition(Math.max(current.minPosition(), next.minPosition()));
                current.setMaxPosition(Math.min(current.maxPosition(), next.maxPosition()));
                positions.remove(i + 1);
            }
            else if(current.chromosome().equals(next.chromosome()))
            {
                current.setMaxPosition(Math.min(current.maxPosition(), next.position()));
                next.setMinPosition(Math.max(next.minPosition(), current.position()));
                i++;
            }
            else
            {
                i++;
            }
        }

        return positions;
    }

}
