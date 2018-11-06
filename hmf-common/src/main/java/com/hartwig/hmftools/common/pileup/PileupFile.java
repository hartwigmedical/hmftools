package com.hartwig.hmftools.common.pileup;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public final class PileupFile {

    private static final String DELIMITER = "\t";

    @NotNull
    public static List<Pileup> read(@NotNull final String filename) throws IOException {
        return Files.readAllLines(new File(filename).toPath()).stream().map(PileupFile::fromString).collect(Collectors.toList());
    }

    @NotNull
    public static Pileup fromString(@NotNull final String line) {

        final Map<String, Integer> insertionMap = Maps.newHashMap();
        final Map<String, Integer> deletionMap = Maps.newHashMap();

        String[] values = line.split(DELIMITER);

        int referenceCount = 0;
        int gCount = 0;
        int aCount = 0;
        int tCount = 0;
        int cCount = 0;
        String indel;
        char prevBase;

        if (values.length >= 5) {
            final String referenceString = values[2];
            final char refBase = referenceString.charAt(0);
            final String readBases = values[4];

            for (int i = 0; i < readBases.length(); i++) {
                switch (Character.toUpperCase(readBases.charAt(i))) {
                    case '+':
                        prevBase = Character.toUpperCase(readBases.charAt(i -1));
                        if (isRef(prevBase)) {
                            prevBase = refBase;
                            referenceCount--;
                        }
                        indel = indel(i, readBases);
                        insertionMap.merge(prevBase + indel, 1, (x, y) -> x + y);
                        i += indel.length() + Integer.toString(indel.length()).length();
                        break;
                    case '-':
                        prevBase = Character.toUpperCase(readBases.charAt(i -1));
                        if (isRef(prevBase)) {
                            prevBase = refBase;
                            referenceCount--;
                        }
                        indel = indel(i, readBases);
                        deletionMap.merge(prevBase + indel, 1, (x, y) -> x + y);
                        i += indel.length() + Integer.toString(indel.length()).length();
                        break;
                    case '^':
                        i++;
                        break;
                    case ',':
                    case '.':
                        referenceCount++;
                        break;
                    case 'G':
                        gCount++;
                        break;
                    case 'A':
                        aCount++;
                        break;
                    case 'T':
                        tCount++;
                        break;
                    case 'C':
                        cCount++;
                        break;
                    default:
                }
            }
        }

        return ImmutablePileup.builder()
                .chromosome(values[0])
                .position(Long.valueOf(values[1]))
                .referenceBase(values[2])
                .readCount(Integer.valueOf(values[3]))
                .referenceCount(referenceCount)
                .gMismatchCount(gCount)
                .aMismatchCount(aCount)
                .tMismatchCount(tCount)
                .cMismatchCount(cCount)
                .insertionCounts(insertionMap)
                .deletionCounts(deletionMap)
                .build();
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

    static boolean inframe(int indelSize) {
        return indelSize % 3 == 0;
    }

    static boolean isRef(char base) {
        return base == ',' || base == '.';
    }
}
