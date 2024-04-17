package com.hartwig.hmftools.sage.old;

import static java.lang.Math.max;

import static com.hartwig.hmftools.sage.old.MicrohomologyContextBuilder.expandMicrohomologyRepeats;
import static com.hartwig.hmftools.sage.old.MicrohomologyContextBuilder.microhomologyAtDeleteFromReadSequence;
import static com.hartwig.hmftools.sage.old.MicrohomologyContextBuilder.microhomologyAtInsert;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;

import java.util.Optional;

import com.hartwig.hmftools.sage.common.RefSequence;

import org.apache.logging.log4j.util.Strings;

import htsjdk.samtools.SAMRecord;

public class ReadContextFactory
{
    private static final int MIN_REPEAT_COUNT = 3;

    private final int mFlankSize;

    public ReadContextFactory(final int flankSize)
    {
        mFlankSize = flankSize;
    }

    public ReadContext createDelContext(
            final String ref, int refPosition, int readIndex, final byte[] readBases, final RefSequence refSequence)
    {
        int refIndex = refSequence.index(refPosition);

        final MicrohomologyContext microhomologyContext = microhomologyAtDeleteFromReadSequence(readIndex, ref, readBases);
        final MicrohomologyContext microhomologyContextWithRepeats = expandMicrohomologyRepeats(microhomologyContext);

        int startIndex = microhomologyContextWithRepeats.position() - MIN_CORE_DISTANCE;
        int length = max(microhomologyContext.length(), microhomologyContextWithRepeats.length() - ref.length() + 1) + 1;
        int endIndex = max(
                microhomologyContextWithRepeats.position() + MIN_CORE_DISTANCE, microhomologyContextWithRepeats.position() + length);

        final Optional<RepeatContext> refRepeatContext = RepeatContextFactory.repeats(refIndex + 1, refSequence.Bases);
        if(refRepeatContext.filter(x -> x.count() >= MIN_REPEAT_COUNT).isPresent())
        {
            final RepeatContext repeat = refRepeatContext.get();
            int repeatStartIndexInReadSpace = repeat.StartIndex - refIndex + readIndex;
            int repeatEndIndexInReadSpace = repeat.EndIndex - refIndex + readIndex;
            startIndex = Math.min(startIndex, repeatStartIndexInReadSpace - 1);
            endIndex = max(endIndex, repeatEndIndexInReadSpace + 1);
        }

        final Optional<RepeatContext> readRepeatContext = RepeatContextFactory.repeats(readIndex + 1, readBases);
        if(readRepeatContext.filter(x -> x.count() >= MIN_REPEAT_COUNT).isPresent())
        {
            final RepeatContext repeat = readRepeatContext.get();
            startIndex = Math.min(startIndex, repeat.StartIndex - 1);
            endIndex = max(endIndex, repeat.EndIndex + 1);
        }

        return ReadContext.fromReadRecord(
                microhomologyContext.toString(),
                readRepeatContext.map(RepeatContext::count).orElse(0),
                readRepeatContext.map(RepeatContext::sequence).orElse(Strings.EMPTY),
                refPosition, readIndex, startIndex, endIndex, mFlankSize, readBases);
    }

    public ReadContext createInsertContext(
            final String alt, int refPosition, int readIndex, final byte[] readBases, final RefSequence refSequence)
    {
        int refIndex = refSequence.index(refPosition);

        final MicrohomologyContext microhomologyContext = microhomologyAtInsert(readIndex, alt.length(), readBases);
        final MicrohomologyContext microhomologyContextWithRepeats = expandMicrohomologyRepeats(microhomologyContext);

        int startIndex = microhomologyContextWithRepeats.position() - MIN_CORE_DISTANCE;
        int length = max(microhomologyContextWithRepeats.length() + 1, alt.length());
        int endIndex = max(
                microhomologyContextWithRepeats.position() + MIN_CORE_DISTANCE, microhomologyContextWithRepeats.position() + length);

        final Optional<RepeatContext> refRepeatContext = RepeatContextFactory.repeats(refIndex + 1, refSequence.Bases);
        if(refRepeatContext.filter(x -> x.count() >= MIN_REPEAT_COUNT).isPresent())
        {
            final RepeatContext repeat = refRepeatContext.get();
            int repeatStartIndexInReadSpace = repeat.StartIndex - refIndex + readIndex;
            int repeatEndIndexInReadSpace = repeat.EndIndex - refIndex + readIndex;
            startIndex = Math.min(startIndex, repeatStartIndexInReadSpace - 1);
            endIndex = max(endIndex, repeatEndIndexInReadSpace + 1);
        }

        final Optional<RepeatContext> readRepeatContext = RepeatContextFactory.repeats(readIndex + 1, readBases);
        if(readRepeatContext.filter(x -> x.count() >= MIN_REPEAT_COUNT).isPresent())
        {
            final RepeatContext repeat = readRepeatContext.get();
            startIndex = Math.min(startIndex, repeat.StartIndex - 1);
            endIndex = max(endIndex, repeat.EndIndex + 1);
        }

        // ensure that MH hasn't reduced the right core index too much
        endIndex = max(endIndex, readIndex + alt.length() - 1 + MIN_CORE_DISTANCE);

        return ReadContext.fromReadRecord(
                microhomologyContext.toString(),
                readRepeatContext.map(RepeatContext::count).orElse(0),
                readRepeatContext.map(RepeatContext::sequence).orElse(Strings.EMPTY),
                refPosition, readIndex, startIndex, endIndex, mFlankSize, readBases);
    }

    public ReadContext createSNVContext(int refPosition, int readIndex, final SAMRecord record, final RefSequence refSequence)
    {
        return createMNVContext(refPosition, readIndex, 1, record.getReadBases(), refSequence);
    }

    public ReadContext createMNVContext(int refPosition, int readIndex, int length, final byte[] readBases, final RefSequence refBases)
    {
        int refIndex = refBases.index(refPosition);
        int startIndex = readIndex - MIN_CORE_DISTANCE;
        int endIndex = readIndex + length - 1 + MIN_CORE_DISTANCE;

        final Optional<RepeatContext> refPriorRepeatContext = RepeatContextFactory.repeats(refIndex - 1, refBases.Bases);
        if(refPriorRepeatContext.filter(x -> x.count() >= MIN_REPEAT_COUNT).isPresent())
        {
            final RepeatContext repeat = refPriorRepeatContext.get();
            int repeatStartIndexInReadSpace = repeat.StartIndex - refIndex + readIndex;
            int repeatEndIndexInReadSpace = repeat.EndIndex - refIndex + readIndex;
            startIndex = Math.min(startIndex, repeatStartIndexInReadSpace - 1);
            endIndex = max(endIndex, repeatEndIndexInReadSpace + 1);
        }

        final Optional<RepeatContext> refPostRepeatContext = RepeatContextFactory.repeats(refIndex + length, refBases.Bases);
        if(refPostRepeatContext.filter(x -> x.count() >= MIN_REPEAT_COUNT).isPresent())
        {
            final RepeatContext repeat = refPostRepeatContext.get();
            int repeatStartIndexInReadSpace = repeat.StartIndex - refIndex + readIndex;
            int repeatEndIndexInReadSpace = repeat.EndIndex - refIndex + readIndex;
            startIndex = Math.min(startIndex, repeatStartIndexInReadSpace - 1);
            endIndex = max(endIndex, repeatEndIndexInReadSpace + 1);
        }

        final Optional<RepeatContext> readRepeatContext = RepeatContextFactory.repeats(readIndex, readBases);
        if(readRepeatContext.filter(x -> x.count() >= MIN_REPEAT_COUNT).isPresent())
        {
            final RepeatContext repeat = readRepeatContext.get();
            startIndex = Math.min(startIndex, repeat.StartIndex - 1);
            endIndex = max(endIndex, repeat.EndIndex + 1);
        }

        return ReadContext.fromReadRecord(
                Strings.EMPTY,
                readRepeatContext.map(RepeatContext::count).orElse(0),
                readRepeatContext.map(RepeatContext::sequence).orElse(Strings.EMPTY),
                refPosition, readIndex, startIndex, endIndex, mFlankSize, readBases);
    }
}
