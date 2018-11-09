package com.hartwig.hmftools.common.pileup;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class PileupFile {

    private static final String DELIMITER = "\t";

    @NotNull
    public static List<Pileup> read(@NotNull final String filename) throws IOException {
        return Files.readAllLines(new File(filename).toPath()).stream().map(PileupFile::fromString).collect(Collectors.toList());
    }

    @NotNull
    public static Pileup fromString(@NotNull final String line) {

        final Map<String, Integer> countMap = Maps.newHashMap();
        final Map<String, Integer> qualityMap = Maps.newHashMap();

        final String[] values = line.split(DELIMITER);
        final int expectedReadCount = Integer.valueOf(values[3]);

        int actualReadCount = 0;
        String indel;
        char prevBase;

        String qualityBases = values.length > 5 ? values[5] : null;

        if (values.length >= 5) {
            final String referenceString = values[2];
            final char refBase = referenceString.charAt(0);
            final String readBases = values[4];

            for (int i = 0; i < readBases.length(); i++) {
                char base = Character.toUpperCase(readBases.charAt(i));
                switch (base) {
                    case '+':
                    case '-':
                        prevBase = Character.toUpperCase(readBases.charAt(i - 1));
                        if (isRef(prevBase)) {
                            prevBase = refBase;
                            countMap.merge(String.valueOf(refBase), -1, (oldValue, value) -> Math.max(0, oldValue + value));
                            qualityMap.merge(String.valueOf(refBase),
                                    -qualityScore(actualReadCount - 1, qualityBases),
                                    (oldValue, value) -> Math.max(0, oldValue + value));
                        }
                        indel = indel(i, readBases);
                        qualityMap.merge(String.valueOf(base) + prevBase + indel,
                                qualityScore(actualReadCount - 1, qualityBases),
                                (oldValue, value) -> oldValue + value);
                        countMap.merge(String.valueOf(base) + prevBase + indel, 1, (oldValue, value) -> oldValue + value);
                        i += indel.length() + Integer.toString(indel.length()).length();
                        break;
                    case ',':
                    case '.':
                        countMap.merge(String.valueOf(refBase), 1, (oldValue, value) -> oldValue + value);
                        qualityMap.merge(String.valueOf(refBase), qualityScore(actualReadCount, qualityBases), (x, y) -> x + y);
                        actualReadCount++;
                        break;
                    case 'G':
                    case 'A':
                    case 'T':
                    case 'C':
                        countMap.merge(String.valueOf(base), 1, (x, y) -> x + y);
                        qualityMap.merge(String.valueOf(base), qualityScore(actualReadCount, qualityBases), (x, y) -> x + y);
                        actualReadCount++;
                        break;
                    case '^':
                        i++;
                        break;
                    case '*':
                        actualReadCount++;
                        break;
                    default:
                }
            }
        }

        if (expectedReadCount != actualReadCount) {
            System.out.println(line);
        }

        return ImmutablePileup.builder()
                .chromosome(values[0])
                .position(Long.valueOf(values[1]))
                .referenceBase(values[2])
                .readCount(expectedReadCount)
                .scoreMap(qualityMap)
                .countMap(countMap)
                .build();
    }

    @VisibleForTesting
    static int qualityScore(int index, @Nullable final String quality) {
        if (quality == null || index >= quality.length()) {
            return 0;
        }

        return (int) quality.charAt(index) - 33;
    }

    @NotNull
    private static String indel(int index, String readBases) {
        int indelSize = indelSize(index + 1, readBases);
        int indelStart = index + 1 + Integer.toString(indelSize).length();
        int indelEnd = indelStart + indelSize;
        return readBases.substring(indelStart, indelEnd).toUpperCase();
    }

    @VisibleForTesting
    static int indelSize(int indelIndex, @NotNull final String text) {
        for (int i = indelIndex + 1; i < text.length(); i++) {
            if (!Character.isDigit(text.charAt(i))) {
                return Integer.valueOf(text.substring(indelIndex, i));
            }
        }

        return Integer.valueOf(text.substring(indelIndex, text.length()));
    }

    private static boolean isRef(char base) {
        return base == ',' || base == '.';
    }
}
