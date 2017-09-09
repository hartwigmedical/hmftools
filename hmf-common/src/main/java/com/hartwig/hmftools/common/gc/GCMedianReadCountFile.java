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
    private static final String EXTENSION = ".purple.gc.median";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static GCMedianReadCount read(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(@NotNull final String fileName, @NotNull final GCMedianReadCount gcMedianReadCount) throws IOException {
        Files.write(new File(fileName).toPath(), toLines(gcMedianReadCount));
    }

    @NotNull
    static GCMedianReadCount fromLines(@NotNull final List<String> lines) throws IOException {

        final int mean = 0;
        final int median = 0;

        final Map<GCBucket, Integer> medianPerBucket = Maps.newHashMap();
        for (int i = 3; i < lines.size(); i++) {
            String[] line = lines.get(i).split(DELIMITER);
            if (line.length == 2) {
                medianPerBucket.put(new ImmutableGCBucket(Integer.valueOf(line[0])), Integer.valueOf(line[1]));
            }
        }

        return new GCMedianReadCountImpl(mean, median, medianPerBucket);
    }

    @NotNull
    static List<String> toLines(@NotNull final GCMedianReadCount gcMedianReadCount) {
        final List<String> lines = Lists.newArrayList();
        lines.add("#SampleMean" + DELIMITER + "SampleMedian");
        lines.add(gcMedianReadCount.meanReadCount() + DELIMITER + gcMedianReadCount.medianReadCount());
        lines.add("#GCBucket" + DELIMITER + "Median");
        for (int i = 0; i <= 100; i++) {
            final GCBucket bucket = new ImmutableGCBucket(i);
            int readCount = gcMedianReadCount.medianReadCount(bucket);
            lines.add(i + "\t" + readCount);
        }
        return lines;
    }
}
