package com.hartwig.hmftools.common.utils.pcf;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
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
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.logging.log4j.util.Strings;

public final class PCFFile
{
    private static final String HEADER_PREFIX = "sampleID";
    private static final String RATIO_EXTENSION = ".cobalt.ratio.pcf";
    private static final String BAF_EXTENSION = ".amber.baf.pcf";

    private static final int COL_CHROMOSOME = 1;
    private static final int COL_POS_START = 3;
    private static final int COL_POS_END = 4;

    public static String generateRatioFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + RATIO_EXTENSION;
    }

    public static String generateBAFFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + BAF_EXTENSION;
    }

    public static ListMultimap<Chromosome, PcfSegment> readPcfFile(final String path)
    {
        ListMultimap<Chromosome, PcfSegment> result = ArrayListMultimap.create();
        try(DelimFileReader dfr = new DelimFileReader(path))
        {
            dfr.stream().forEach(row ->
            {
                final String chrName = row.get(0);
                Chromosome chromosome = HumanChromosome.fromString(chrName);
                int start = row.getInt(1);
                int end = row.getInt(2);
                double meanRatio = row.getDouble(3);
                result.put(chromosome, new PcfSegment(chrName, start, end, meanRatio));
            });
        }
        return result;
    }

    public static void write(final String filename, final GenomeIntervals data)
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

    public static Map<Chromosome,List<PCFPosition>> readPositions(int windowSize, PCFSource source, String filename) throws IOException
    {
        Map<Chromosome,List<PCFPosition>> chrPositionsMap = Maps.newHashMap();
        Window window = new Window(windowSize);

        PCFPosition lastPcfPosition = null;

        String currentChromosome = Strings.EMPTY;
        int minPosition = 1;
        List<PCFPosition> pcfPositions = null;

        List<String> lines = Files.readAllLines(new File(filename).toPath());
        String header = lines.get(0);

        boolean inOldFormat = header.startsWith(HEADER_PREFIX);
        int chrIndex = inOldFormat ? COL_CHROMOSOME : 0;
        int posStartIndex = inOldFormat ? COL_POS_START : 1;
        int posEndIndex = inOldFormat ? COL_POS_END : 2;

        for(int i = 1; i < lines.size(); i++)
        {
            String line = lines.get(i);
            String[] values = line.split(TSV_DELIM);

            String chrStr = values[chrIndex];

            if(!HumanChromosome.contains(chrStr))
                continue;

            // each entry results in 2 PCF positions since gaps are converted into entries as well
            if(!chrStr.equals(currentChromosome))
            {
                currentChromosome = chrStr;
                pcfPositions = Lists.newArrayList();
                chrPositionsMap.put(HumanChromosome.fromString(chrStr), pcfPositions);

                minPosition = 1; // reset for start of new chromosome
                lastPcfPosition = null;
            }

            int rawStart = Integer.parseInt(values[posStartIndex]);
            int rawEnd = Integer.parseInt(values[posEndIndex]);
            int start = inOldFormat ? window.start(rawStart) : rawStart;
            int end = inOldFormat ? window.start(rawEnd) + windowSize : rawEnd + 1;

            if(lastPcfPosition != null)
            {
                lastPcfPosition.setMaxPosition(start); // mark the boundary for the previous PCF position
            }

            PCFPosition pcfStart = new PCFPosition(source, chrStr, start);
            pcfStart.setMinPosition(minPosition); // taken from the previous position's end
            pcfStart.setMaxPosition(start);
            pcfPositions.add(pcfStart);

            minPosition = end;

            PCFPosition pcfEnd = new PCFPosition(source, chrStr, end);
            pcfEnd.setMinPosition(end);
            pcfEnd.setMaxPosition(end);
            pcfPositions.add(pcfEnd);

            lastPcfPosition = pcfEnd;
        }

        for(List<PCFPosition> positions : chrPositionsMap.values())
        {
            mergePositions(positions);
        }

        return chrPositionsMap;
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
