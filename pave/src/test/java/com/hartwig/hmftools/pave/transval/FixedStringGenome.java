package com.hartwig.hmftools.pave.transval;

import org.jetbrains.annotations.NotNull;

public class FixedStringGenome extends SimpleTestGenome
{
    private final String Bases;
    public FixedStringGenome(@NotNull final String bases)
    {
        this.Bases = bases;
    }

    @Override
    public String getBaseString(final String chromosome, final int posStart, final int posEnd)
    {
        return Bases.substring(posStart, posEnd);
    }

    @Override
    public byte[] getBases(final String chromosome, final int posStart, final int posEnd)
    {
        return getBaseString(chromosome, posStart, posEnd).getBytes();
    }
}
