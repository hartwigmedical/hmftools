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
    public static ReadContext createDelContext(@NotNull final String ref, int refPosition, int readIndex, @NotNull final SAMRecord record,
            int refIndex, byte[] refBases) {

        final MicrohomologyContext microhomologyContext = microhomologyAtDeleteFromReadSequence(readIndex, ref, record.getReadBases());
        final MicrohomologyContext microhomologyContextWithRepeats = expandMicrohomologyRepeats(microhomologyContext);

        int startIndex = microhomologyContextWithRepeats.position();
        int length = Math.max(microhomologyContext.length(), microhomologyContextWithRepeats.length() - ref.length() + 1) + 1;
        int endIndex = Math.min(record.getReadBases().length, startIndex + length);

        final Optional<RepeatContext> repeatContext = RepeatContextFactory.repeats(readIndex + 1, record.getReadBases());
        if (repeatContext.isPresent()) {
            final RepeatContext repeat = repeatContext.get();
            startIndex = Math.min(startIndex, repeat.startIndex() - 1);
            endIndex = Math.max(endIndex, repeat.endIndex() + 1);
        }

        return new ReadContext(microhomologyContext.toString(),
                repeatContext.map(RepeatContext::count).orElse(0),
                repeatContext.map(RepeatContext::sequence).orElse(Strings.EMPTY),
                refPosition,
                readIndex,
                Math.max(startIndex, 0),
                Math.min(endIndex, record.getReadBases().length - 1),
                DEFAULT_BUFFER,
                refBases,
                record);
    }

    @NotNull
    public static ReadContext createInsertContext(@NotNull final String alt, int refPosition, int readIndex,
            @NotNull final SAMRecord record, final int refIndex, byte[] refBases) {

        final MicrohomologyContext microhomologyContext = microhomologyAtInsert(readIndex, alt.length(), record.getReadBases());
        final MicrohomologyContext microhomologyContextWithRepeats = expandMicrohomologyRepeats(microhomologyContext);

        int startIndex = microhomologyContextWithRepeats.position();
        int length = Math.max(microhomologyContextWithRepeats.length() + 1, alt.length());
        int endIndex = Math.min(record.getReadBases().length, startIndex + length);

        final Optional<RepeatContext> repeatContext = RepeatContextFactory.repeats(readIndex + 1, record.getReadBases());
        if (repeatContext.isPresent()) {
            final RepeatContext repeat = repeatContext.get();
            startIndex = Math.min(startIndex, repeat.startIndex() - 1);
            endIndex = Math.max(endIndex, repeat.endIndex() + 1);
        }

        return new ReadContext(microhomologyContext.toString(),
                repeatContext.map(RepeatContext::count).orElse(0),
                repeatContext.map(RepeatContext::sequence).orElse(Strings.EMPTY),
                refPosition,
                readIndex,
                Math.max(startIndex, 0),
                Math.min(endIndex, record.getReadBases().length - 1),
                DEFAULT_BUFFER,
                refBases,
                record);
    }

    @NotNull
    public static ReadContext createSNVContext(int refPosition, int readIndex, @NotNull final SAMRecord record, int refIndex,
            byte[] refBases) {

        int startIndex = readIndex;
        int endIndex = readIndex;
        final Optional<RepeatContext> repeatContext = RepeatContextFactory.repeats(readIndex, record.getReadBases());
        if (repeatContext.isPresent()) {
            final RepeatContext repeat = repeatContext.get();
            startIndex = Math.min(startIndex, repeat.startIndex() - 1);
            endIndex = Math.max(endIndex, repeat.endIndex() + 1);
        }

        return new ReadContext(Strings.EMPTY,
                repeatContext.map(RepeatContext::count).orElse(0),
                repeatContext.map(RepeatContext::sequence).orElse(Strings.EMPTY),
                refPosition,
                readIndex,
                Math.max(startIndex, 0),
                Math.min(endIndex, record.getReadBases().length - 1),
                DEFAULT_BUFFER,
                refBases,
                record);
    }

    @NotNull
    public static ReadContext dummy(int refPosition, @NotNull final String alt) {
        return new ReadContext(Strings.EMPTY, refPosition, 0, 0, 0, DEFAULT_BUFFER, alt.getBytes());
    }

}
