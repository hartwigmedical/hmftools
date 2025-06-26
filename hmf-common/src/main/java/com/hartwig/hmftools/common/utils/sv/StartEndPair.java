package com.hartwig.hmftools.common.utils.sv;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.Objects;

public class StartEndPair<A>
{
    public final A Start;
    public final A End;

    public StartEndPair(A start, A end)
    {
        Start = start;
        End = end;
    }

    public A get(int seIndex) { return seIndex == SE_START ? Start : End; }
    public A get(boolean isStart) { return isStart ? Start : End; }
    public A start() { return Start; }
    public A end() { return End; }

    public String toString()
    {
        return "Pair[" + Start + "," + End + "]";
    }

    @Override
    public boolean equals(Object other)
    {
        return other instanceof StartEndPair<?>
                && Objects.equals(Start, ((StartEndPair)other).Start)
                && Objects.equals(End, ((StartEndPair)other).End);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Start, End);
    }

    public static <A,B> StartEndPair<A> of(A a, A b)
    {
        return new StartEndPair<>(a,b);
    }
}
