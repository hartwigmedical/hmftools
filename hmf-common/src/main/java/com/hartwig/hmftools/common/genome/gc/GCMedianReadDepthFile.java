package com.hartwig.hmftools.common.genome.gc;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public final class GCMedianReadDepthFile
{
    private static final String EXTENSION = ".cobalt.gc.median.tsv";
    private static final int ASSUMED_READ_LENGTH = 151;

    @NotNull
    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static GCMedianReadDepth read(final String filename) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(final String fileName, final GCMedianReadDepth gcMedianReadDepth) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(gcMedianReadDepth));
    }

    @NotNull
    private static GCMedianReadDepth fromLines(final List<String> lines)
    {
        boolean useReadDepth = true;
        double mean = 0;
        double median = 0;
        int i = 0;

        // read the first line, see if it uses upper case #SampleMean or lower case #sampleMean
        // lower case sampleMean is newer version that uses read depth
        // upper case SampleMean is old version that uses read count
        if(!lines.isEmpty())
        {
             if(lines.get(i).startsWith("#SampleMean"))
             {
                useReadDepth = false;
             }
        }

        // skip the #sampleMean     sampleMedian line
        ++i;

        if(lines.size() > i)
        {
            String[] line = lines.get(i++).split(TSV_DELIM);
            mean = Double.parseDouble(line[0]);
            median = Double.parseDouble(line[1]);
        }

        // skip the #gcBucket       median line
        ++i;

        Map<GCBucket, Double> medianPerBucket = new HashMap<>();

        for(; i < lines.size(); i++)
        {
            String[] line = lines.get(i).split(TSV_DELIM);
            if(line.length == 2)
            {
                medianPerBucket.put(new ImmutableGCBucket(Integer.parseInt(line[0])), Double.parseDouble(line[1]));
            }
        }

        if(!useReadDepth)
        {
            // if we are parsing old version of cobalt output, convert it all to depth by assuming read length of 151
            mean = convertReadCount(mean);
            median = convertReadCount(median);
            medianPerBucket = medianPerBucket.entrySet().stream().collect(Collectors.toMap(
                    Map.Entry::getKey,
                    entry -> convertReadCount(entry.getValue())));
        }

        return new GCMedianReadDepth(mean, median, medianPerBucket);
    }

    @NotNull
    private static List<String> toLines(final GCMedianReadDepth gcMedianReadDepth)
    {
        final List<String> lines = new ArrayList<>();
        lines.add("#sampleMean" + TSV_DELIM + "sampleMedian");
        lines.add(String.format("%.2f" + TSV_DELIM + "%.2f", gcMedianReadDepth.meanReadDepth(), gcMedianReadDepth.medianReadDepth()));
        lines.add("#gcBucket" + TSV_DELIM + "median");
        for(int i = 0; i <= 100; i++)
        {
            final GCBucket bucket = new ImmutableGCBucket(i);
            double readDepth = gcMedianReadDepth.medianReadDepth(bucket);
            if(readDepth > 0)
            {
                lines.add(String.format("%d" + TSV_DELIM + "%.2f", i, readDepth));
            }
        }
        return lines;
    }

    // this is backward compatibility conversion from read count to read depth
    // It assumes read length is 151.
    private static double convertReadCount(final double readCount)
    {
        if(readCount <= 0)
            return readCount;

        return readCount * ASSUMED_READ_LENGTH / GCProfileFactory.WINDOW_SIZE;
    }
}
