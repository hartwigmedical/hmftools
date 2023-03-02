package com.hartwig.hmftools.datamodel.genome.refgenome;
import org.jetbrains.annotations.NotNull;

public enum RefGenomeVersion
{
    V37("37", true),
    V38("38", false);

    @NotNull
    private final String mIdentifier;
    private final boolean mIs37;

    RefGenomeVersion(@NotNull final String identifier, final boolean is37)
    {
        mIdentifier = identifier;
        mIs37 = is37;
    }

    public boolean is37() { return mIs37; }
    public boolean is38 () { return !mIs37; }

    public String identifier() { return mIdentifier; }
}
