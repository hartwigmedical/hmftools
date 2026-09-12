package com.hartwig.hmftools.common.cobalt;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.gc.GCBucket;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCBucket;

public final class CobaltGcMedianFile
{
    private static final String EXTENSION = ".cobalt.gc.median.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static GcMedianReadDepth read(final String filename) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(final String fileName, final GcMedianReadDepth gcMedianReadDepth) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(gcMedianReadDepth));
    }

    private static GcMedianReadDepth fromLines(final List<String> lines)
    {
        double mean = 0;
        double median = 0;
        int i = 0;

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

        return new GcMedianReadDepth(mean, median, medianPerBucket);
    }

    private static List<String> toLines(final GcMedianReadDepth gcMedianReadDepth)
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
}
