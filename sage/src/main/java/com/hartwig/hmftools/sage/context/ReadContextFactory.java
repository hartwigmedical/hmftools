package com.hartwig.hmftools.sage.context;

import static com.hartwig.hmftools.common.variant.Microhomology.expandMicrohomologyRepeats;
import static com.hartwig.hmftools.common.variant.Microhomology.microhomologyAtDelete;
import static com.hartwig.hmftools.common.variant.Microhomology.microhomologyAtInsert;

import com.hartwig.hmftools.common.variant.MicrohomologyContext;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContextFactory {

    private static final int DEFAULT_BUFFER = 25;

    @NotNull
    public static ReadContextImproved createDelContext(@NotNull final String ref, int refPosition, int readIndex,
            @NotNull final SAMRecord record, int refIndex, byte[] refBases) {

        final MicrohomologyContext microhomologyContext = microhomologyAtDelete(refIndex, ref.length(), refBases);
        final MicrohomologyContext microhomologyContextWithRepeats = expandMicrohomologyRepeats(microhomologyContext);

        int startIndex = microhomologyContextWithRepeats.position();
        int length = Math.max(microhomologyContext.length(), microhomologyContextWithRepeats.length() - ref.length() + 1) + 1;
        int endIndex = Math.min(record.getReadBases().length, startIndex + length);

        return new ReadContextImproved(0, refPosition, readIndex, startIndex, endIndex, DEFAULT_BUFFER, refBases, record);
    }

    @NotNull
    public static ReadContextImproved createInsertContext(@NotNull final String alt, int refPosition, int readIndex,
            @NotNull final SAMRecord record, byte[] refBases) {

        final MicrohomologyContext microhomologyContext = microhomologyAtInsert(readIndex, alt.length(), record.getReadBases());
        final MicrohomologyContext microhomologyContextWithRepeats = expandMicrohomologyRepeats(microhomologyContext);

        int startIndex = microhomologyContextWithRepeats.position();
        int length = Math.max(microhomologyContextWithRepeats.length() + 1, alt.length());

        int endIndex = Math.min(record.getReadBases().length, startIndex + length);

        return new ReadContextImproved(0, refPosition, readIndex, startIndex, endIndex, DEFAULT_BUFFER, refBases, record);
    }

    @NotNull
    public static ReadContextImproved createSNVContext(int refPosition, int readIndex, @NotNull final SAMRecord record, byte[] refBases) {
        return new ReadContextImproved(0, refPosition, readIndex, readIndex, readIndex, DEFAULT_BUFFER, refBases, record);
    }

    @NotNull
    public static ReadContextImproved dummy(int refPosition, @NotNull final String alt) {
        return new ReadContextImproved(0, refPosition, 0, 0, 0, DEFAULT_BUFFER, alt.getBytes());
    }

}
