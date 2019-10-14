package com.hartwig.hmftools.sage.context;

import com.hartwig.hmftools.common.variant.Microhomology;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContextFactory {

    private static final int DEFAULT_BUFFER = 25;

    @NotNull
    public static ReadContextImproved createDelContext(@NotNull final String ref, int refPosition, int readIndex,
            @NotNull final SAMRecord record, int refIndex, byte[] refBases) {
        final String microhomology = Microhomology.microhomologyAtDelete(refIndex, new String(refBases), ref);
        return new ReadContextImproved(refPosition,
                readIndex,
                readIndex - 1,
                readIndex + 1 + microhomology.length(),
                DEFAULT_BUFFER,
                record.getReadBases());
    }

    @NotNull
    public static ReadContextImproved createInsertContext(@NotNull final String alt, int refPosition, int readIndex,
            @NotNull final SAMRecord record, int refIndex, byte[] refBases) {
        final String microhomology = Microhomology.microhomologyAtInsert(refIndex, new String(refBases), alt);
        return new ReadContextImproved(refPosition,
                readIndex,
                readIndex - 1,
                readIndex + 1 + microhomology.length() + alt.length() - 1,
                DEFAULT_BUFFER,
                record.getReadBases());
    }

    @NotNull
    public static ReadContextImproved createSNVContext(int refPosition, int readIndex, @NotNull final SAMRecord record) {

        return createOld(refPosition, readIndex, record);
        //        int leftIndex = readIndex;
        //        int rightIndex = readIndex;
        //
        //        byte[] readBases = record.getReadBases();
        //        byte readIndexBase = readBases[readIndex];
        //        for (int i = 1; i <= readIndex; i++) {
        //            if (readBases[leftIndex - 1] != readIndexBase) {
        //                break;
        //            }
        //            leftIndex--;
        //        }
        //
        //        for (int i = readIndex + 1; i < readBases.length; i++) {
        //            if (readBases[i] != readIndexBase) {
        //                break;
        //            }
        //            rightIndex++;
        //        }
        //
        //        return new ReadContextImproved(refPosition, readIndex, leftIndex, rightIndex, DEFAULT_BUFFER, readBases);
    }

    @NotNull
    public static ReadContextImproved createOld(int refPosition, int readIndex, @NotNull final SAMRecord record) {
        return new ReadContextImproved(refPosition, readIndex, readIndex, readIndex, DEFAULT_BUFFER, record.getReadBases());
    }

    @NotNull
    public static ReadContextImproved dummy(int refPosition, @NotNull final String alt) {
        return new ReadContextImproved(refPosition, 0, 0, 0, DEFAULT_BUFFER, alt.getBytes());
    }

}
