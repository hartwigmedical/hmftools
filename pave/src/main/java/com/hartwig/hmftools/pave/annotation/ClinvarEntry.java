package com.hartwig.hmftools.pave.annotation;

public class ClinvarEntry
{
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final String Significance;
    public final String Conflict;

    public ClinvarEntry(final int position, final String ref, final String alt, final String significance, final String conflict)
    {
        Position = position;
        Ref = ref;
        Alt = alt;
        Significance = significance;
        Conflict = conflict;
    }

    public boolean matches(final int position, final String ref, final String alt)
    {
        return Position == position && Ref.equals(ref) && Alt.equals(alt);
    }

    public String toString()
    {
        return String.format("%d %s>%s details(%s - %s)",
                Position, Ref, Alt, Significance, Conflict);
    }
}
