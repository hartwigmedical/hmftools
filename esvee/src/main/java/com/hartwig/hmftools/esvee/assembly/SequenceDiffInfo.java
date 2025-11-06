package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.BASE;

import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;

public class SequenceDiffInfo
{
    public final int ReadIndex; // index within the read's bases
    public final int ConsensusIndex; // index of this read's mismatch in the assembly, ie not in the read itself

    public final String Bases; // either the alt SNV base, the novel indel or the repeat sequence if a contraction or expansion
    public final SequenceDiffType Type;

    // repeat-diff fields
    public final int RepeatCount; // repeat count in the read
    public final int RepeatIndexBegin; // index where the repeat begins in the read, from the build direction

    public final BaseQualType QualType;
    public double MismatchPenalty;

    public static final SequenceDiffInfo UNSET = new SequenceDiffInfo(
            -1, -1, "", SequenceDiffType.UNSET, BaseQualType.LOW, 0, -1);

    public SequenceDiffInfo(
            int readIndex, int consensusIndex, final String bases, final SequenceDiffType type, final BaseQualType qualType)
    {
        this(readIndex, consensusIndex, bases, type, qualType, 0, -1);
    }

    public SequenceDiffInfo(
            int readIndex, int consensusIndex, final String bases, final SequenceDiffType type, final BaseQualType qualType,
            int repeatCount, int repeatIndexBegin)
    {
        ReadIndex = readIndex;
        ConsensusIndex = consensusIndex;
        Bases = bases;
        Type = type;
        QualType = qualType;
        RepeatCount = repeatCount;
        RepeatIndexBegin = repeatIndexBegin;

        MismatchPenalty = 0;
    }

    public static SequenceDiffInfo fromSnv(final ReadParseState read, int consensusIndex)
    {
        return new SequenceDiffInfo(
                read.readIndex(), consensusIndex, String.valueOf((char)read.currentBase()), BASE, read.qualType(read.readIndex()));
    }

    public static SequenceDiffInfo fromMatch(int readIndex, int consensusIndex)
    {
        return new SequenceDiffInfo(readIndex, consensusIndex, null, SequenceDiffType.MATCH, BaseQualType.HIGH);
    }

    // repeat related methods
    public int repeatIndex(final RepeatInfo consensusRepeat, final boolean buildForwards, boolean getStart)
    {
        if(buildForwards == getStart)
            return RepeatIndexBegin;

        if(buildForwards)
            return RepeatIndexBegin + RepeatCount * consensusRepeat.repeatLength() - 1;
        else
            return RepeatIndexBegin - RepeatCount * consensusRepeat.repeatLength() + 1;
    }

    public static int repeatIndex(
            final int readIndex, final int readRepeatCount, final RepeatInfo consensusRepeat, final int previousRepeatBases,
            final boolean buildForwards, boolean getStart)
    {
        if(buildForwards == getStart)
        {
            return buildForwards ? readIndex - previousRepeatBases : readIndex + previousRepeatBases;
        }

        if(buildForwards)
        {
            int repeatIndexStart = readIndex - previousRepeatBases;
            return repeatIndexStart + readRepeatCount * consensusRepeat.repeatLength() - 1;
        }
        else
        {
            int repeatIndexEnd = readIndex + previousRepeatBases;
            return repeatIndexEnd - readRepeatCount * consensusRepeat.repeatLength() + 1;
        }
    }

    public String toString()
    {
        if(Type == SequenceDiffType.MATCH)
            return format("index(r=%d c=%d): %s", ReadIndex, ConsensusIndex, Type);

        String info = format("index(r=%d c=%d): %s %s penalty(%.2f)", ReadIndex, ConsensusIndex, Type, Bases, MismatchPenalty);
        return RepeatCount == 0 ? info : format("%s x%d)", info, RepeatCount);
    }
}
