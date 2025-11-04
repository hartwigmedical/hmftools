package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.BASE;

public class SequenceDiffInfo
{
    public final int ReadIndex; // index within the read's bases
    public final int ConsensusIndex; // index of this read's mismatch in the assembly, ie not in the read itself

    public final String Bases; // either the alt SNV base, the novel indel or the repeat sequence if a contraction or expansion
    public final SequenceDiffType Type;
    public final int RepeatCount; // only applicable for repeats, repeat count in the read

    public final BaseQualType QualType;
    public double MismatchPenalty;

    public static final SequenceDiffInfo UNSET = new SequenceDiffInfo(
            -1, -1, "", SequenceDiffType.UNSET, 0, BaseQualType.LOW);

    public SequenceDiffInfo(
            int readIndex, int consensusIndex, final String bases, final SequenceDiffType type, int repeatCount, final BaseQualType qualType)
    {
        ReadIndex = readIndex;
        ConsensusIndex = consensusIndex;
        Bases = bases;
        Type = type;
        RepeatCount = repeatCount;
        QualType = qualType;
        MismatchPenalty = 0;
    }

    public static SequenceDiffInfo fromSnv(final ReadParseState read, int consensusIndex)
    {
        return new SequenceDiffInfo(
                read.readIndex(), consensusIndex, String.valueOf((char)read.currentBase()), BASE, 0,
                read.qualType(read.readIndex()));
    }

    public static SequenceDiffInfo fromMatch(int readIndex, int consensusIndex)
    {
        return new SequenceDiffInfo(readIndex, consensusIndex, null, SequenceDiffType.MATCH, (short)0, BaseQualType.HIGH);
    }

    public String toString()
    {
        if(Type == SequenceDiffType.MATCH)
            return format("index(r=%d c=%d): %s", ReadIndex, ConsensusIndex, Type);

        String info = format("index(r=%d c=%d): %s %s penalty(%g)", ReadIndex, ConsensusIndex, Type, Bases, MismatchPenalty);
        return RepeatCount == 0 ? info : format("%s x%d)", info, RepeatCount);
    }
}
