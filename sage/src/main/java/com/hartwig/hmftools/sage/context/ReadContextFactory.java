package com.hartwig.hmftools.sage.context;

import static com.hartwig.hmftools.common.variant.Microhomology.expandMicrohomologyRepeats;
import static com.hartwig.hmftools.common.variant.Microhomology.microhomologyAtDeleteFromReadSequence;
import static com.hartwig.hmftools.common.variant.Microhomology.microhomologyAtInsert;

import java.util.Optional;

import com.hartwig.hmftools.common.variant.MicrohomologyContext;
import com.hartwig.hmftools.common.variant.repeat.RepeatContext;
import com.hartwig.hmftools.common.variant.repeat.RepeatContextFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContextFactory {

    private static final int DEFAULT_BUFFER = 25;

    @NotNull
    public static ReadContextImproved createDelContext(@NotNull final String ref, int refPosition, int readIndex,
            @NotNull final SAMRecord record, int refIndex, byte[] refBases) {

        final String refString = new String(refBases);
        final String readString = new String(record.getReadBases());

        final MicrohomologyContext microhomologyContext = microhomologyAtDeleteFromReadSequence(readIndex, ref, record.getReadBases());
        final MicrohomologyContext microhomologyContextWithRepeats = expandMicrohomologyRepeats(microhomologyContext);

        int startIndex = microhomologyContextWithRepeats.position();
        int length = Math.max(microhomologyContext.length(), microhomologyContextWithRepeats.length() - ref.length() + 1) + 1;
        int endIndex = Math.min(record.getReadBases().length, startIndex + length);

        final Optional<RepeatContext> repeatContext = RepeatContextFactory.repeats(readIndex + 1, readString);
        if (repeatContext.isPresent()) {
            final RepeatContext repeat = repeatContext.get();
            startIndex = Math.min(startIndex, readIndex - repeat.backwardsCount() * repeat.sequence().length());
            endIndex = Math.max(endIndex,
                    Math.min(record.getReadBases().length - 1, readIndex + repeat.forwardsCount() * repeat.sequence().length() + 1));
        }

        final String repeat = repeatContext.map(RepeatContext::sequence).orElse(Strings.EMPTY);
        return new ReadContextImproved(microhomologyContext.toString(),
                repeat,
                refPosition,
                readIndex,
                startIndex,
                endIndex,
                DEFAULT_BUFFER,
                refBases,
                record);
    }

    @NotNull
    public static ReadContextImproved createInsertContext(@NotNull final String alt, int refPosition, int readIndex,
            @NotNull final SAMRecord record, final int refIndex, byte[] refBases) {

        final String refString = new String(refBases);
        final String readString = new String(record.getReadBases());

        final MicrohomologyContext microhomologyContext = microhomologyAtInsert(readIndex, alt.length(), record.getReadBases());
        final MicrohomologyContext microhomologyContextWithRepeats = expandMicrohomologyRepeats(microhomologyContext);

        int startIndex = microhomologyContextWithRepeats.position();
        int length = Math.max(microhomologyContextWithRepeats.length() + 1, alt.length());
        int endIndex = Math.min(record.getReadBases().length, startIndex + length);


        final Optional<RepeatContext> repeatContext = RepeatContextFactory.repeats(readIndex + 1, readString);
        if (repeatContext.isPresent()) {
            final RepeatContext repeat = repeatContext.get();
            startIndex = Math.min(startIndex, readIndex - repeat.backwardsCount() * repeat.sequence().length());
            endIndex = Math.max(endIndex,
                    Math.min(record.getReadBases().length - 1, readIndex + repeat.forwardsCount() * repeat.sequence().length() + 1));
        }

        final String repeat = repeatContext.map(RepeatContext::sequence).orElse(Strings.EMPTY);
        return new ReadContextImproved(microhomologyContext.toString(),
                repeat,
                refPosition,
                readIndex,
                startIndex,
                endIndex,
                DEFAULT_BUFFER,
                refBases,
                record);
    }

    @NotNull
    public static ReadContextImproved createSNVContext(int refPosition, int readIndex, @NotNull final SAMRecord record, int refIndex,
            byte[] refBases) {
        return new ReadContextImproved(Strings.EMPTY,
                Strings.EMPTY,
                refPosition,
                readIndex,
                readIndex,
                readIndex,
                DEFAULT_BUFFER,
                refBases,
                record);
    }

    @NotNull
    public static ReadContextImproved dummy(int refPosition, @NotNull final String alt) {
        return new ReadContextImproved(Strings.EMPTY, Strings.EMPTY, refPosition, 0, 0, 0, DEFAULT_BUFFER, alt.getBytes());
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
