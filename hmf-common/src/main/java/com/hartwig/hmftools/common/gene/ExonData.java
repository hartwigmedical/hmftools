package com.hartwig.hmftools.common.gene;

public class ExonData
{
    public final int TransId;
    public final int Start;
    public final int End;
    public final int Rank;
    public final int PhaseStart;
    public final int PhaseEnd;

    public ExonData(int transId, int start, int end, int rank, int phaseStart, int phaseEnd)
    {
        TransId = transId;
        Start = start;
        End = end;
        Rank = rank;
        PhaseStart = phaseStart;
        PhaseEnd = phaseEnd;
    }

    public int length() { return End - Start; }
    public int baseLength() { return length() + 1; }

    public String toString()
    {
        return String.format("%d: range(%d -> %d) phases(%d -> %d)", Rank, Start, End, PhaseStart, PhaseEnd);
    }

}
