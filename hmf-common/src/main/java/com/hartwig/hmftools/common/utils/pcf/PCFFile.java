package com.hartwig.hmftools.common.utils.pcf;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.genome.region.Window;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PCFFile
{
    private static final String HEADER_PREFIX = "sampleID";
    private static final String RATIO_EXTENSION = ".cobalt.ratio.pcf";
    private static final String COBALT_PCF_EXTENSION = ".cobalt.pcf.tsv";
    private static final String BAF_EXTENSION = ".amber.baf.pcf";

    private static final int COL_CHROMOSOME = 1;
    private static final int COL_POS_START = 3;
    private static final int COL_POS_END = 4;

    @NotNull
    public static String generateCobaltPcfFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + COBALT_PCF_EXTENSION;
    }

    @NotNull
    public static String generateRatioFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + RATIO_EXTENSION;
    }

    @NotNull
    public static String generateBAFFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + BAF_EXTENSION;
    }

    public static ListMultimap<Chromosome, ChrBaseRegion> readCobaltPcfFile(String path)
    {
        ListMultimap<Chromosome, ChrBaseRegion> result = ArrayListMultimap.create();
        try(DelimFileReader dfr = new DelimFileReader(path))
        {
            dfr.stream().forEach(row ->
            {
                final String chrName = row.get(0);
                Chromosome chromosome = HumanChromosome.fromString(chrName);
                int start = row.getInt(1);
                int end = row.getInt(2);
                result.put(chromosome, new ChrBaseRegion(chrName, start, end));
            });
        }
        return result;
    }

    public static void write(String filename, GenomeIntervals data)
    {
        List<ChrBaseRegion> ratios = data.regionsList();
        List<String> columns = List.of("sampleID", "chrom", "arm", "start.pos", "end.pos", "n.probes", "mean");
        DelimFileWriter.write(filename, columns, ratios,
                (ratio, row) ->
                {
                    row.set(columns.get(0), "unused");
                    row.set(columns.get(1), ratio.chromosome());
                    row.set(columns.get(2), "unused");
                    row.set(columns.get(3), ratio.start());
                    row.set(columns.get(4), ratio.end());
                    row.set(columns.get(5), -1);
                    row.set(columns.get(6), -1);
                });
    }

    @NotNull
    public static ListMultimap<Chromosome, PCFPosition> readPositions(int windowSize, final PCFSource source, final String filename)
            throws IOException
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
                String[] values = line.split(TSV_DELIM);
                String chromosomeName = values[COL_CHROMOSOME];
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

                    int start = window.start(Integer.parseInt(values[COL_POS_START]));
                    int end = window.start(Integer.parseInt(values[COL_POS_END])) + windowSize;
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

    public static Map<String, List<BaseRegion>> loadChrBaseRegions(final String filename)
    {
        if(filename == null)
        {
            return Collections.emptyMap();
        }

        try
        {
            Map<String, List<BaseRegion>> regionsMap = Maps.newHashMap();

            List<BaseRegion> regions = null;
            String currentChromosome = "";

            List<String> lines = Files.readAllLines(new File(filename).toPath());
            lines.remove(0);

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                String chromosome = values[COL_CHROMOSOME];
                int posStart = Integer.parseInt(values[COL_POS_START]);
                int posEnd = Integer.parseInt(values[COL_POS_END]);

                if(!currentChromosome.equals(chromosome))
                {
                    currentChromosome = chromosome;
                    regions = regionsMap.get(chromosome);

                    if(regions == null)
                    {
                        regions = Lists.newArrayList();
                        regionsMap.put(chromosome, regions);
                    }
                }

                regions.add(new BaseRegion(posStart, posEnd));
            }

            return regionsMap;
        }
        catch(IOException e)
        {
            return null;
        }
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
        String[] values = line.split(TSV_DELIM);
        String chromosome = values[COL_CHROMOSOME];
        int posStart = Integer.parseInt(values[COL_POS_START]);
        int posEnd = Integer.parseInt(values[COL_POS_END]) + windowSize - 1;
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
