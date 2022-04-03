package com.hartwig.hmftools.common.genome.gc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class GCMedianReadCountFile {

    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".cobalt.gc.median.tsv";

    private GCMedianReadCountFile() {
    }

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
        double mean = 0;
        double median = 0;
        if (lines.size() >= 2) {
            String[] line = lines.get(1).split(DELIMITER);
            mean = Double.parseDouble(line[0]);
            median = Double.parseDouble(line[1]);
        }

        final Map<GCBucket, Double> medianPerBucket = new HashMap<>();
        if (extendRange) {
            if (lines.size() > 3) {
                String[] minLine = lines.get(3).split(DELIMITER);
                double minMedian = Double.parseDouble(minLine[1]);
                for (int i = 0; i < Double.parseDouble(minLine[0]); i++) {
                    medianPerBucket.put(new ImmutableGCBucket(i), minMedian);
                }

                String[] maxLine = lines.get(lines.size() - 1).split(DELIMITER);
                double maxMedian = Double.parseDouble(maxLine[1]);
                for (int i = Integer.parseInt(maxLine[0]) + 1; i <= 100; i++) {
                    medianPerBucket.put(new ImmutableGCBucket(i), maxMedian);
                }
            }
        }

        for (int i = 3; i < lines.size(); i++) {
            String[] line = lines.get(i).split(DELIMITER);
            if (line.length == 2) {
                medianPerBucket.put(new ImmutableGCBucket(Integer.parseInt(line[0])), Double.parseDouble(line[1]));
            }
        }

        return new GCMedianReadCount(mean, median, medianPerBucket);
    }

    @NotNull
    private static List<String> toLines(@NotNull final GCMedianReadCount gcMedianReadCount) {
        final List<String> lines = Lists.newArrayList();
        lines.add("#SampleMean" + DELIMITER + "SampleMedian");
        lines.add(String.format("%.2f" + DELIMITER + "%.2f", gcMedianReadCount.meanReadCount(), gcMedianReadCount.medianReadCount()));
        lines.add("#GCBucket" + DELIMITER + "Median");
        for (int i = 0; i <= 100; i++) {
            final GCBucket bucket = new ImmutableGCBucket(i);
            double readCount = gcMedianReadCount.medianReadCount(bucket);
            if (readCount > 0) {
                lines.add(String.format("%d" + DELIMITER + "%.2f", i, readCount));
            }
        }
        return lines;
    }
}
