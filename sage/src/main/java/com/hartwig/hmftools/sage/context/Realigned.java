package com.hartwig.hmftools.sage.context;

import static com.hartwig.hmftools.sage.context.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.context.RealignedType.SHORTENED;

import java.util.Optional;

import com.hartwig.hmftools.common.variant.repeat.RepeatContext;
import com.hartwig.hmftools.common.variant.repeat.RepeatContextFactory;
import com.hartwig.hmftools.sage.read.ReadContext;

import org.jetbrains.annotations.NotNull;

public class Realigned {

    public static final int MAX_REPEAT_SIZE = 5;
    private static final RealignedContext NONE = new RealignedContext(RealignedType.NONE, 0);
    private static final RealignedContext EXACT = new RealignedContext(RealignedType.EXACT, 0);

    @NotNull
    public static RealignedContext realignedAroundIndex(@NotNull final ReadContext readContext, final int otherBaseIndex, final byte[] otherBases,
            int maxSize) {

        if (realigned(readContext.readBasesPositionIndex(),
                readContext.readBasesLeftFlankIndex(),
                readContext.readBasesRightFlankIndex(),
                readContext.readBases(),
                otherBaseIndex,
                otherBases,
                maxSize)) {
            return EXACT;
        }

        return jitter(readContext.readBasesPositionIndex(),
                readContext.readBasesLeftCentreIndex(),
                readContext.readBasesRightCentreIndex(),
                readContext.readBases(),
                otherBaseIndex,
                otherBases);

    }

    static boolean realigned(int baseIndex, int baseStartIndex, int baseEndIndex, final byte[] bases, final int otherBaseIndex,
            final byte[] otherBases, int maxDistance) {

        for (int offset = -maxDistance; offset <= maxDistance; offset++) {
            if (realigned(baseIndex, baseStartIndex, baseEndIndex, bases, otherBaseIndex + offset, otherBases)) {
                return true;
            }
        }

        return false;
    }

    private static boolean realigned(int baseIndex, int baseStartIndex, int baseEndIndex, final byte[] bases, int otherIndex, byte[] otherBases) {

        int leftOffset = baseIndex - baseStartIndex;
        int otherStartIndex = otherIndex - leftOffset;
        int expectedLength = baseEndIndex - baseStartIndex + 1;
        return otherStartIndex >= 0 && otherStartIndex + expectedLength <= otherBases.length
                && matchingBasesFromLeft(baseStartIndex, baseEndIndex, bases, otherStartIndex, otherBases) == expectedLength;

    }

    @NotNull
    static RealignedContext jitter(int baseIndex, int baseStartIndex, int baseEndIndex, final byte[] bases, int otherIndex,
            byte[] otherBases) {

        int leftOffset = baseIndex - baseStartIndex;
        int otherStartIndex = otherIndex - leftOffset;
        return jitter(baseStartIndex, baseEndIndex, bases, otherStartIndex, otherBases);

    }

    @NotNull
    private static RealignedContext jitter(int baseStartIndex, int baseEndIndex, final byte[] bases, int otherStartIndex,
            byte[] otherBases) {
        int exactLength = baseEndIndex - baseStartIndex + 1;

        int matchingBases = matchingBasesFromLeft(baseStartIndex, baseEndIndex, bases, otherStartIndex, otherBases);
        if (matchingBases == exactLength) {
            return EXACT;
        }

        int baseNextIndex = baseStartIndex + matchingBases;
        int otherNextIndex = otherStartIndex + matchingBases;

        final Optional<RepeatContext> optionalRepeat = RepeatContextFactory.repeats(otherNextIndex - 1, otherBases);
        if (!optionalRepeat.isPresent()) {
            return NONE;
        }

        final RepeatContext repeat = optionalRepeat.get();

        int repeatLength = repeat.sequence().length();
        if (repeatLength == 0) {
            return NONE;
        }

        int matchingBasesLengthened = matchingBasesFromLeft(baseNextIndex - repeatLength, baseEndIndex, bases, otherNextIndex, otherBases);
        if (matchingBasesLengthened > 0 && matchingBases + matchingBasesLengthened == exactLength + repeatLength) {
            return new RealignedContext(LENGTHENED, repeat.count());
        }

        int matchingBasesShortened = matchingBasesFromLeft(baseNextIndex + repeatLength, baseEndIndex, bases, otherNextIndex, otherBases);
        if (matchingBasesShortened > 0 && matchingBases + matchingBasesShortened == exactLength - repeatLength) {
            return new RealignedContext(SHORTENED, repeat.count());
        }

        return NONE;
    }

    private static int matchingBasesFromLeft(int startIndex, int endIndex, byte[] bases, int otherIndex, byte[] otherBases) {
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

}
