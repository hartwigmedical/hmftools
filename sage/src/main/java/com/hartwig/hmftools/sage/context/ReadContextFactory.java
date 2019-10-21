package com.hartwig.hmftools.sage.context;

import static com.hartwig.hmftools.common.variant.Microhomology.expandMicrohomologyRepeats;
import static com.hartwig.hmftools.common.variant.Microhomology.microhomologyAtDelete;
import static com.hartwig.hmftools.common.variant.Microhomology.microhomologyAtInsert;

import java.util.Optional;

import com.hartwig.hmftools.common.variant.MicrohomologyContext;
import com.hartwig.hmftools.common.variant.repeat.RepeatContext;
import com.hartwig.hmftools.common.variant.repeat.RepeatContextFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContextFactory {

    private static final int REPEAT_BASES = 10;
    private static final int DEFAULT_BUFFER = 25;

    @NotNull
    public static ReadContextImproved createDelContext(@NotNull final String ref, int refPosition, int readIndex,
            @NotNull final SAMRecord record, int refIndex, byte[] refBases) {

        final MicrohomologyContext microhomologyContext = microhomologyAtDelete(refIndex, ref.length(), refBases);
        final MicrohomologyContext microhomologyContextWithRepeats = expandMicrohomologyRepeats(microhomologyContext);

        int startIndex = microhomologyContextWithRepeats.position();
        int length = Math.max(microhomologyContext.length(), microhomologyContextWithRepeats.length() - ref.length() + 1) + 1;
        int endIndex = Math.min(record.getReadBases().length, startIndex + length);

        final Optional<RepeatContext> repeatContext =
                RepeatContextFactory.repeats(0, new String(refBases, refIndex + 1, Math.min(REPEAT_BASES, refBases.length - refIndex - 1)));
        int jitter = repeatContext.map(x -> jitter(microhomologyContext, x)).orElse(0);

        return new ReadContextImproved(jitter, refPosition, readIndex, startIndex, endIndex, DEFAULT_BUFFER, refBases, record);
    }

    @NotNull
    public static ReadContextImproved createInsertContext(@NotNull final String alt, int refPosition, int readIndex,
            @NotNull final SAMRecord record, final int refIndex, byte[] refBases) {

        final MicrohomologyContext microhomologyContext = microhomologyAtInsert(readIndex, alt.length(), record.getReadBases());
        final MicrohomologyContext microhomologyContextWithRepeats = expandMicrohomologyRepeats(microhomologyContext);

        int startIndex = microhomologyContextWithRepeats.position();
        int length = Math.max(microhomologyContextWithRepeats.length() + 1, alt.length());

        int endIndex = Math.min(record.getReadBases().length, startIndex + length);

        final Optional<RepeatContext> repeatContext =
                RepeatContextFactory.repeats(0, new String(refBases, refIndex + 1, Math.min(REPEAT_BASES, refBases.length - refIndex - 1)));
        int jitter = repeatContext.map(x -> jitter(microhomologyContext, x)).orElse(0);

        return new ReadContextImproved(jitter, refPosition, readIndex, startIndex, endIndex, DEFAULT_BUFFER, refBases, record);
    }

    @NotNull
    public static ReadContextImproved createSNVContext(int refPosition, int readIndex, @NotNull final SAMRecord record, byte[] refBases) {
        return new ReadContextImproved(0, refPosition, readIndex, readIndex, readIndex, DEFAULT_BUFFER, refBases, record);
    }

    @NotNull
    public static ReadContextImproved dummy(int refPosition, @NotNull final String alt) {
        return new ReadContextImproved(0, refPosition, 0, 0, 0, DEFAULT_BUFFER, alt.getBytes());
    }

    private static int jitter(@NotNull final MicrohomologyContext microhomologyContext, @NotNull final RepeatContext context) {
        if (context.sequence().length() > microhomologyContext.length()) {
            return 0;
        }

        byte[] repeatBytes = context.sequence().getBytes();
        for (int i = 0; i < repeatBytes.length; i++) {
            if (repeatBytes[repeatBytes.length - 1 - i] != microhomologyContext.readSequence()[
                    microhomologyContext.homologyIndex() + microhomologyContext.length() - 1 - i]) {
                return 0;
            }
        }

        return context.sequence().length();
    }

}
