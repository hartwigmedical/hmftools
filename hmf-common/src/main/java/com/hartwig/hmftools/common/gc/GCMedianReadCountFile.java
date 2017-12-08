package com.hartwig.hmftools.common.gc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public enum GCMedianReadCountFile {
    ;

    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".cobalt.gc.median";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static GCMedianReadCount read(boolean extendRange, @NotNull final String filename) throws IOException {
        return fromLines(extendRange, Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(@NotNull final String fileName, @NotNull final GCMedianReadCount gcMedianReadCount) throws IOException {
        Files.write(new File(fileName).toPath(), toLines(gcMedianReadCount));
    }

    @NotNull
    private static GCMedianReadCount fromLines(boolean extendRange, @NotNull final List<String> lines) {
        int mean = 0;
        int median = 0;
        if (lines.size() >= 2) {
            String[] line = lines.get(1).split(DELIMITER);
            mean = Integer.valueOf(line[0]);
            median = Integer.valueOf(line[1]);
        }

        final Map<GCBucket, Integer> medianPerBucket = Maps.newHashMap();
        if (extendRange) {
            if (lines.size() > 3) {
                String[] minLine = lines.get(3).split(DELIMITER);
                int minMedian = Integer.valueOf(minLine[1]);
                for (int i = 0; i < Integer.valueOf(minLine[0]); i++) {
                    medianPerBucket.put(new ImmutableGCBucket(i), minMedian);
                }

                String[] maxLine = lines.get(lines.size() - 1).split(DELIMITER);
                int maxMedian = Integer.valueOf(maxLine[1]);
                for (int i = Integer.valueOf(maxLine[0]) + 1; i <= 100; i++) {
                    medianPerBucket.put(new ImmutableGCBucket(i), maxMedian);
                }
            }
        }

        for (int i = 3; i < lines.size(); i++) {
            String[] line = lines.get(i).split(DELIMITER);
            if (line.length == 2) {
                medianPerBucket.put(new ImmutableGCBucket(Integer.valueOf(line[0])), Integer.valueOf(line[1]));
            }
        }

        return new GCMedianReadCountImpl(mean, median, medianPerBucket);
    }

    @NotNull
    private static List<String> toLines(@NotNull final GCMedianReadCount gcMedianReadCount) {
        final List<String> lines = Lists.newArrayList();
        lines.add("#SampleMean" + DELIMITER + "SampleMedian");
        lines.add(gcMedianReadCount.meanReadCount() + DELIMITER + gcMedianReadCount.medianReadCount());
        lines.add("#GCBucket" + DELIMITER + "Median");
        for (int i = 0; i <= 100; i++) {
            final GCBucket bucket = new ImmutableGCBucket(i);
            int readCount = gcMedianReadCount.medianReadCount(bucket);
            if (readCount > 0) {
                lines.add(i + "\t" + readCount);
            }
        }
        return lines;
    }
}
