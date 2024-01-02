package com.hartwig.hmftools.sage.candidate_;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.variant.Microhomology.expandMicrohomologyRepeats;
import static com.hartwig.hmftools.common.variant.Microhomology.microhomologyAtDeleteFromReadSequence;
import static com.hartwig.hmftools.common.variant.Microhomology.microhomologyAtInsert;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;

import java.util.Optional;

import com.hartwig.hmftools.common.variant.MicrohomologyContext;
import com.hartwig.hmftools.common.variant.repeat.RepeatContext;
import com.hartwig.hmftools.common.variant.repeat.RepeatContextFactory_;

import org.apache.logging.log4j.util.Strings;

import htsjdk.samtools.SAMRecord;

public class ReadContextFactory_
{
    private static final int MIN_REPEAT_COUNT = 3;

    private final int mFlankSize;

    public ReadContextFactory_(final int flankSize)
    {
        mFlankSize = flankSize;
    }

    /**
     * Creates a read context for an MNV.
     * <p>
     * startIndex and endIndex for read context initially contains the MVN + a buffer of MIN_CORE_DISTANCE.
     * Expands these to contain (+1 base buffer) any repeats, of count at least MIN_REPEAT_COUNT, in the ref before and after the MVN, and
     * also any repeats in the read that start at or before the first MVN base.
     * <p>
     * No homology for an MVN.
     */
    public ReadContext_ createMNVContext(int refPosition, int readIndex, int length, final byte[] readBases, final IndexedBases_ refBases)
    {
        // Index of refBases with position mapped to refPosition.
        int refIndex = refBases.index(refPosition);

        // Pads out by MIN_CORE_DISTANCE.
        int startIndex = readIndex - MIN_CORE_DISTANCE;
        int endIndex = readIndex + length - 1 + MIN_CORE_DISTANCE;

        // Get repeats starting before refIndex in refBases.
        // If there is a repeat with repeats at least MIN_REPEAT_COUNT, we contain the repeat.
        final Optional<RepeatContext> refPriorRepeatContext = RepeatContextFactory_.repeats(refIndex - 1, refBases.Bases);
        if(refPriorRepeatContext.filter(x -> x.count() >= MIN_REPEAT_COUNT).isPresent())
        {
            final RepeatContext repeat = refPriorRepeatContext.get();
            int repeatStartIndexInReadSpace = repeat.startIndex() - refIndex + readIndex;
            int repeatEndIndexInReadSpace = repeat.endIndex() - refIndex + readIndex;
            startIndex = Math.min(startIndex, repeatStartIndexInReadSpace - 1);
            endIndex = max(endIndex, repeatEndIndexInReadSpace + 1);
        }

        // Gets repeats starting after the MVN in refBases.
        // If there is a repeat with repeats at least MIN_REPEAT_COUNT, we contain the repeat + 1 base for buffer.
        final Optional<RepeatContext> refPostRepeatContext = RepeatContextFactory_.repeats(refIndex + length, refBases.Bases);
        if(refPostRepeatContext.filter(x -> x.count() >= MIN_REPEAT_COUNT).isPresent())
        {
            final RepeatContext repeat = refPostRepeatContext.get();
            int repeatStartIndexInReadSpace = repeat.startIndex() - refIndex + readIndex;
            int repeatEndIndexInReadSpace = repeat.endIndex() - refIndex + readIndex;
            startIndex = Math.min(startIndex, repeatStartIndexInReadSpace - 1);
            endIndex = max(endIndex, repeatEndIndexInReadSpace + 1);
        }

        // Get repeats starting at readIndex in the read (first base of MNV).
        // If there is a repeat with repeats at least MIN_REPEAT_COUNT, we contain the repeat + 1 base for buffer.
        final Optional<RepeatContext> readRepeatContext = RepeatContextFactory_.repeats(readIndex, readBases);
        if(readRepeatContext.filter(x -> x.count() >= MIN_REPEAT_COUNT).isPresent())
        {
            final RepeatContext repeat = readRepeatContext.get();
            startIndex = Math.min(startIndex, repeat.startIndex() - 1);
            endIndex = max(endIndex, repeat.endIndex() + 1);
        }

        // No microhomology.
        return ReadContext_.fromReadRecord(
                Strings.EMPTY,
                readRepeatContext.map(RepeatContext::count).orElse(0),
                readRepeatContext.map(RepeatContext::sequence).orElse(Strings.EMPTY),
                refPosition, readIndex, startIndex, endIndex, mFlankSize, readBases);
    }

    public ReadContext_ createSNVContext(int refPosition, int readIndex, final SAMRecord record, final IndexedBases_ refBases)
    {
        return createMNVContext(refPosition, readIndex, 1, record.getReadBases(), refBases);
    }

    // TODO
    /**
     * Creates a read context for a DEL.
     * <p>
     * startIndex and endIndex for read context initially contains the microhomology repeats and a buffer of MIN_CORE_DISTANCE.
     * Expands these to contain (+1 base buffer) any repeats, of count at least MIN_REPEAT_COUNT, in the ref starting at or before the first
     * deleted base, and also any repeats in the read that start at or before first base at or after the DEL.
     */
    public ReadContext_ createDelContext(
            final String ref, int refPosition, int readIndex, final byte[] readBases, final IndexedBases_ refBases)
    {
        // Index of refBases with position mapped to refPosition.
        int refIndex = refBases.index(refPosition);

        // TODO: Computes the microhomology.
        final MicrohomologyContext microhomologyContext = microhomologyAtDeleteFromReadSequence(readIndex, ref, readBases);
        final MicrohomologyContext microhomologyContextWithRepeats = expandMicrohomologyRepeats(microhomologyContext);

        // TODO: The initial startIndex and endIndex contains the microhomology repeats and a buffer of MIN_CORE_DISTANCE.
        int startIndex = microhomologyContextWithRepeats.position() - MIN_CORE_DISTANCE;
        int length = max(microhomologyContext.length(), microhomologyContextWithRepeats.length() - ref.length() + 1) + 1;
        int endIndex = max(
                microhomologyContextWithRepeats.position() + MIN_CORE_DISTANCE, microhomologyContextWithRepeats.position() + length);

        // Gets repeats starting at or before the first deleted base.
        // If there is a repeat with repeats at least MIN_REPEAT_COUNT, we contain the repeat + 1 base for buffer.
        final Optional<RepeatContext> refRepeatContext = RepeatContextFactory_.repeats(refIndex + 1, refBases.Bases);
        if(refRepeatContext.filter(x -> x.count() >= MIN_REPEAT_COUNT).isPresent())
        {
            final RepeatContext repeat = refRepeatContext.get();
            int repeatStartIndexInReadSpace = repeat.startIndex() - refIndex + readIndex;
            int repeatEndIndexInReadSpace = repeat.endIndex() - refIndex + readIndex;
            startIndex = Math.min(startIndex, repeatStartIndexInReadSpace - 1);
            endIndex = max(endIndex, repeatEndIndexInReadSpace + 1);
        }

        // Get repeats starting at or before readIndex + 1 in the read (first base after the DEL).
        // If there is a repeat with repeats at least MIN_REPEAT_COUNT, we contain the repeat + 1 base for buffer.
        final Optional<RepeatContext> readRepeatContext = RepeatContextFactory_.repeats(readIndex + 1, readBases);
        if(readRepeatContext.filter(x -> x.count() >= MIN_REPEAT_COUNT).isPresent())
        {
            final RepeatContext repeat = readRepeatContext.get();
            startIndex = Math.min(startIndex, repeat.startIndex() - 1);
            endIndex = max(endIndex, repeat.endIndex() + 1);
        }

        return ReadContext_.fromReadRecord(
                microhomologyContext.toString(),
                readRepeatContext.map(RepeatContext::count).orElse(0),
                readRepeatContext.map(RepeatContext::sequence).orElse(Strings.EMPTY),
                refPosition, readIndex, startIndex, endIndex, mFlankSize, readBases);
    }

    // TODO
    /**
     * Creates a read context for an INS.
     * <p>
     * startIndex and endIndex for read context initially contains the microhomology repeats and a buffer of MIN_CORE_DISTANCE.
     * Expands these to contain (+1 base buffer) any repeats, of count at least MIN_REPEAT_COUNT, in the ref starting at or before the first
     * deleted base, and also any repeats in the read that start at or before first base at or after the DEL.
     * <p>
     * Finally exapnds endIndex to contain the full insert + MIN_CORE_DISTANCE.
     */
    public ReadContext_ createInsertContext(
            final String alt, int refPosition, int readIndex, final byte[] readBases, final IndexedBases_ refBases)
    {
        // Index of refBases with position mapped to refPosition.
        int refIndex = refBases.index(refPosition);

        // TODO: Computes the microhomology.
        final MicrohomologyContext microhomologyContext = microhomologyAtInsert(readIndex, alt.length(), readBases);
        final MicrohomologyContext microhomologyContextWithRepeats = expandMicrohomologyRepeats(microhomologyContext);

        // TODO: The initial startIndex and endIndex contains the microhomology repeats and a buffer of MIN_CORE_DISTANCE.
        int startIndex = microhomologyContextWithRepeats.position() - MIN_CORE_DISTANCE;
        int length = max(microhomologyContextWithRepeats.length() + 1, alt.length());
        int endIndex = max(
                microhomologyContextWithRepeats.position() + MIN_CORE_DISTANCE, microhomologyContextWithRepeats.position() + length);

        // Gets repeats starting at or before the first deleted base.
        // If there is a repeat with repeats at least MIN_REPEAT_COUNT, we contain the repeat + 1 base for buffer.
        final Optional<RepeatContext> refRepeatContext = RepeatContextFactory_.repeats(refIndex + 1, refBases.Bases);
        if(refRepeatContext.filter(x -> x.count() >= MIN_REPEAT_COUNT).isPresent())
        {
            final RepeatContext repeat = refRepeatContext.get();
            int repeatStartIndexInReadSpace = repeat.startIndex() - refIndex + readIndex;
            int repeatEndIndexInReadSpace = repeat.endIndex() - refIndex + readIndex;
            startIndex = Math.min(startIndex, repeatStartIndexInReadSpace - 1);
            endIndex = max(endIndex, repeatEndIndexInReadSpace + 1);
        }

        // Get repeats starting at or before readIndex + 1 in the read (first base after the DEL).
        // If there is a repeat with repeats at least MIN_REPEAT_COUNT, we contain the repeat + 1 base for buffer.
        final Optional<RepeatContext> readRepeatContext = RepeatContextFactory_.repeats(readIndex + 1, readBases);
        if(readRepeatContext.filter(x -> x.count() >= MIN_REPEAT_COUNT).isPresent())
        {
            final RepeatContext repeat = readRepeatContext.get();
            startIndex = Math.min(startIndex, repeat.startIndex() - 1);
            endIndex = max(endIndex, repeat.endIndex() + 1);
        }

        // ensure that MH hasn't reduced the right core index too much
        endIndex = max(endIndex, readIndex + alt.length() - 1 + MIN_CORE_DISTANCE);

        return ReadContext_.fromReadRecord(
                microhomologyContext.toString(),
                readRepeatContext.map(RepeatContext::count).orElse(0),
                readRepeatContext.map(RepeatContext::sequence).orElse(Strings.EMPTY),
                refPosition, readIndex, startIndex, endIndex, mFlankSize, readBases);
    }
}
