package com.hartwig.hmftools.sage.context;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContextFactory {

    private static final int DEFAULT_BUFFER = 25;

    @NotNull
    public static ReadContextImproved create(int refPosition, int readIndex, @NotNull final SAMRecord record) {

        int leftIndex = readIndex;
        int rightIndex = readIndex;

        byte[] readBases = record.getReadBases();
        byte readIndexBase = readBases[readIndex];
        for (int i = 1; i <= readIndex; i++) {
            if (readBases[leftIndex - 1] != readIndexBase) {
                break;
            }
            leftIndex--;
        }

        for (int i = readIndex + 1; i < readBases.length; i++) {
            if (readBases[i] != readIndexBase) {
                break;
            }
            rightIndex++;
        }

        return new ReadContextImproved(refPosition, leftIndex, rightIndex, DEFAULT_BUFFER, readBases);
    }

    @NotNull
    public static ReadContextImproved createOld(int refPosition, int readIndex, @NotNull final SAMRecord record) {
        return new ReadContextImproved(refPosition, readIndex, readIndex, DEFAULT_BUFFER, record.getReadBases());
    }

}
