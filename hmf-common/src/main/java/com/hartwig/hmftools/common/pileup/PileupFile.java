package com.hartwig.hmftools.common.pileup;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

public class PileupFile {

    private static final String DELIMITER = "\t";

    @NotNull
    public static List<Pileup> read(@NotNull final String filename) throws IOException {
        return Files.readAllLines(new File(filename).toPath()).stream().map(PileupFile::fromString).collect(Collectors.toList());
    }

    @VisibleForTesting
    @NotNull
    static Pileup fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);

        int referenceCount = 0;
        int gCount = 0;
        int aCount = 0;
        int tCount = 0;
        int cCount = 0;
        int insertions = 0;
        int deletions = 0;

        if (values.length >= 5) {
            final String readBases = values[4];
            for (int i = 0; i < readBases.length(); i++) {
                switch (Character.toUpperCase(readBases.charAt(i))) {
                    case '+':
                        insertions++;
                        i += indelStringSize(i + 1, readBases);
                        break;
                    case '-':
                        deletions++;
                        i += indelStringSize(i + 1, readBases);
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
                .insertions(insertions)
                .deletions(deletions)
                .build();
    }

    @VisibleForTesting
    static int indelStringSize(int indelIndex, @NotNull final String text) {
        for (int i = indelIndex + 1; i < text.length(); i++) {
            if (!Character.isDigit(text.charAt(i))) {
                return i - indelIndex + Integer.valueOf(text.substring(indelIndex, i));
            }
        }

        return Integer.valueOf(text.substring(indelIndex, text.length()));
    }
}
