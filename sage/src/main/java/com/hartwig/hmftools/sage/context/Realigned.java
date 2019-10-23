package com.hartwig.hmftools.sage.context;

import static com.hartwig.hmftools.sage.context.RealignedType.EXACT;
import static com.hartwig.hmftools.sage.context.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.context.RealignedType.NONE;
import static com.hartwig.hmftools.sage.context.RealignedType.SHORTENED;

import com.hartwig.hmftools.common.variant.repeat.RepeatContextFactory;

import org.jetbrains.annotations.NotNull;

public class Realigned {

    private static final int MIN_REPEAT_COUNT = 4;
    private static final int MAX_REPEAT_SIZE = 5;

    @NotNull
    public RealignedType realigned(@NotNull final ReadContext readContext, final byte[] otherBases) {
        return realigned(readContext.leftFlankStartIndex(), readContext.rightFlankEndIndex(), readContext.readBases(), otherBases);
    }

    @NotNull
    public RealignedType realigned(int baseStartIndex, int baseEndIndex, final byte[] bases, final byte[] otherBases) {

        int exactLength = baseEndIndex - baseStartIndex + 1;

        for (int i = 0; i <= otherBases.length - exactLength + MAX_REPEAT_SIZE; i++) {
            RealignedType type = realigned(baseStartIndex, baseEndIndex, bases, i, otherBases);
            if (type != NONE) {
                return type;
            }
        }

        return NONE;
    }

    @NotNull
    public RealignedType realigned(int baseStartIndex, int baseEndIndex, final byte[] bases, int otherIndex, byte[] otherBases) {
        int exactLength = baseEndIndex - baseStartIndex + 1;

        int matchingBases = matchingBasesFromLeft(baseStartIndex, baseEndIndex, bases, otherIndex, otherBases);
        if (matchingBases == exactLength) {
            return EXACT;
        }

        if (matchingBases < MIN_REPEAT_COUNT) {
            return NONE;
        }

        int baseNextIndex = baseStartIndex + matchingBases;
        int otherNextIndex = otherIndex + matchingBases;

        int sizeOfRepeatPriorTo = repeatCount(otherNextIndex, otherBases);
        if (sizeOfRepeatPriorTo == 0) {
            return NONE;
        }

        int matchingBasesShortened =
                matchingBasesFromLeft(baseNextIndex + sizeOfRepeatPriorTo, baseEndIndex, bases, otherNextIndex, otherBases);
        if (matchingBasesShortened > 0 && matchingBases + matchingBasesShortened == exactLength - sizeOfRepeatPriorTo) {
            return SHORTENED;
        }

        int matchingBasesLengthened =
                matchingBasesFromLeft(baseNextIndex - sizeOfRepeatPriorTo, baseEndIndex, bases, otherNextIndex, otherBases);
        if (matchingBasesLengthened > 0 && matchingBases + matchingBasesLengthened == exactLength + sizeOfRepeatPriorTo) {
            return LENGTHENED;
        }

        return NONE;
    }

    private int matchingBasesFromLeft(int startIndex, int endIndex, byte[] bases, int otherIndex, byte[] otherBases) {
        if (startIndex < 0) {
            return 0;
        }

        int maxLength = Math.min(endIndex - startIndex + 1, otherBases.length - otherIndex);

        for (int i = 0; i < maxLength; i++) {
            if (bases[startIndex + i] != otherBases[otherIndex + i]) {
                return i;
            }
        }

        return maxLength;
    }

    private int repeatCount(int index, byte[] bases) {
        for (int i = 1; i <= MAX_REPEAT_SIZE; i++) {
            int repeats = RepeatContextFactory.backwardRepeats(index - i, i, bases);
            if (repeats >= MIN_REPEAT_COUNT - 1) {
                return i;
            }

        }

        return 0;

    }

}
