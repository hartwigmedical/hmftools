package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

public class SequenceDiffInfo
{
    public final int ReadIndex; // index within the read's bases
    public final int ConsensusIndex; // index of this read's mismatch in the assembly, ie not in the read itself

    public final String Bases; // either the alt SNV base, the novel indel or the repeat sequence if a contraction or expansion
    public final SequenceDiffType Type;
    public final int RepeatCount; // only applicable for repeats, repeat count in the read

    public double MismatchPenalty;

    public static final SequenceDiffInfo UNSET = new SequenceDiffInfo(
            -1, -1, "", SequenceDiffType.UNSET, 0);

    public SequenceDiffInfo(int readIndex, int consensusIndex, final String bases, final SequenceDiffType type, int repeatCount)
    {
        ReadIndex = readIndex;
        ConsensusIndex = consensusIndex;
        Bases = bases;
        Type = type;
        RepeatCount = repeatCount;
        MismatchPenalty = 0;
    }

    public static SequenceDiffInfo fromSnv(int readIndex, int consensusIndex, final byte base)
    {
        return new SequenceDiffInfo(readIndex, consensusIndex, String.valueOf((char)base), SequenceDiffType.BASE, 0);
    }

    public static SequenceDiffInfo fromMatch(int readIndex, int consensusIndex)
    {
        return new SequenceDiffInfo(readIndex, consensusIndex, null, SequenceDiffType.MATCH, (short)0);
    }

    public String toString()
    {
        if(Type == SequenceDiffType.MATCH)
            return format("index(r=%d c=%d): %s", ReadIndex, ConsensusIndex, Type);

        String info = format("index(r=%d c=%d): %s %s penalty(%g)", ReadIndex, ConsensusIndex, Type, Bases, MismatchPenalty);
        return RepeatCount == 0 ? info : format("%s x%d)", info, RepeatCount);
    }
}
