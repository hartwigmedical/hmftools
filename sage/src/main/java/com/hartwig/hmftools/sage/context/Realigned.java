package com.hartwig.hmftools.sage.context;

import static com.hartwig.hmftools.sage.context.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.context.RealignedType.SHORTENED;

import com.hartwig.hmftools.common.variant.repeat.RepeatContextFactory;
import com.hartwig.hmftools.sage.read.ReadContext;

import org.jetbrains.annotations.NotNull;

public class Realigned {

    private static final int MIN_REPEAT_COUNT = 4;
    private static final int MAX_REPEAT_SIZE = 5;
    private static final RealignedContext NONE = new RealignedContext(RealignedType.NONE, 0);
    private static final RealignedContext EXACT = new RealignedContext(RealignedType.EXACT, 0);
    private static final Repeat NO_REPEAT = new Repeat(0, 0);

    @NotNull
    public RealignedContext realignedInEntireRecord(@NotNull final ReadContext readContext, final byte[] otherBases) {
        return realigned(readContext.leftFlankStartIndex(), readContext.rightFlankEndIndex(), readContext.readBases(), otherBases);
    }

    @NotNull
    public RealignedContext realignedAroundIndex(@NotNull final ReadContext readContext, final int otherBaseIndex, final byte[] otherBases) {
        return realigned(readContext.leftFlankStartIndex(), readContext.rightFlankEndIndex(), readContext.readBases(), otherBaseIndex, otherBases, MAX_REPEAT_SIZE);
    }

    @NotNull
    private RealignedContext realigned(int baseStartIndex, int baseEndIndex, final byte[] bases,  final int otherBaseIndex, final byte[] otherBases, int maxDistance) {

        RealignedContext context = realigned(baseStartIndex, baseEndIndex, bases, otherBaseIndex, otherBases);
        if (context.type() != RealignedType.NONE) {
            return context;
        }

        int exactLength = baseEndIndex - baseStartIndex + 1;

        for (int i = otherBaseIndex + 1; i <= Math.min(otherBaseIndex + maxDistance, otherBases.length - exactLength + MAX_REPEAT_SIZE); i++) {
            context = realigned(baseStartIndex, baseEndIndex, bases, i, otherBases);
            if (context.type() != RealignedType.NONE) {
                return context;
            }
        }

        for (int i = otherBaseIndex - 1; i >= Math.max(otherBaseIndex - maxDistance, 0); i--) {
            context = realigned(baseStartIndex, baseEndIndex, bases, i, otherBases);
            if (context.type() != RealignedType.NONE) {
                return context;
            }
        }

        return NONE;
    }

    @NotNull
    RealignedContext realigned(int baseStartIndex, int baseEndIndex, final byte[] bases, final byte[] otherBases) {

        int exactLength = baseEndIndex - baseStartIndex + 1;

        for (int i = 0; i <= otherBases.length - exactLength + MAX_REPEAT_SIZE; i++) {
            RealignedContext context = realigned(baseStartIndex, baseEndIndex, bases, i, otherBases);
            if (context.type() != RealignedType.NONE) {
                return context;
            }
        }

        return NONE;
    }

    @NotNull
    private RealignedContext realigned(int baseStartIndex, int baseEndIndex, final byte[] bases, int otherIndex, byte[] otherBases) {
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

        final Repeat repeat = repeatCount(otherNextIndex, otherBases);
        int repeatLength = repeat.repeatLength;
        if (repeatLength == 0) {
            return NONE;
        }

        int matchingBasesShortened = matchingBasesFromLeft(baseNextIndex + repeatLength, baseEndIndex, bases, otherNextIndex, otherBases);
        if (matchingBasesShortened > 0 && matchingBases + matchingBasesShortened == exactLength - repeatLength) {
            return new RealignedContext(SHORTENED, repeat.repeatCount);
        }

        int matchingBasesLengthened = matchingBasesFromLeft(baseNextIndex - repeatLength, baseEndIndex, bases, otherNextIndex, otherBases);
        if (matchingBasesLengthened > 0 && matchingBases + matchingBasesLengthened == exactLength + repeatLength) {
            return new RealignedContext(LENGTHENED, repeat.repeatCount + 1);
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

    private Repeat repeatCount(int index, byte[] bases) {
        for (int i = 1; i <= MAX_REPEAT_SIZE; i++) {
            int repeats = RepeatContextFactory.backwardRepeats(index - i, i, bases) + 1;
            if (repeats >= MIN_REPEAT_COUNT) {
                return new Repeat(i, repeats);
            }
        }

        return NO_REPEAT;

    }

    private static class Repeat {
        private final int repeatLength;
        private final int repeatCount;

        Repeat(final int repeatLength, final int repeatCount) {
            this.repeatLength = repeatLength;
            this.repeatCount = repeatCount;
        }
    }

}
