package com.hartwig.hmftools.telo.breakend;

import java.util.Comparator;

import org.jetbrains.annotations.NotNull;

// an added telomere usually means a section of the chromosome is
// broken off
// there might be a case for working out if the telomere is attached to the centremere
// or somewhere else??
// there are two types of such locations
// 1. On the sen
public class TelomericBreakEnd implements Comparable<TelomericBreakEnd>
{
    public enum Type
    {
        // potential new telomere
        RIGHT_G_TELOMERIC,
        LEFT_C_TELOMERIC,
        // potential join with telomere
        RIGHT_C_TELOMERIC,
        LEFT_G_TELOMERIC
    }

    // +1 orientation right side been lost and got a C rich telomere
    // -1 orientation left side been lost and got a G rich telomere
    private final Type mType;
    private final String mChromosome;
    private final int mPosition;

    public TelomericBreakEnd(Type type, String chromosome, int position)
    {
        mType = type;
        mChromosome = chromosome;
        mPosition = position;
    }

    public Type getType() { return mType; }
    public String getChromosome() { return mChromosome; }
    public int getPosition() { return mPosition; }

    // Overriding compareTo() method to allow us to sort
    @Override
    public int compareTo(@NotNull TelomericBreakEnd o)
    {
        return Comparator.comparing(TelomericBreakEnd::getType)
                .thenComparing(TelomericBreakEnd::getChromosome)
                .thenComparingInt(TelomericBreakEnd::getPosition)
                .compare(this, o);
    }

    @Override
    public String toString()
    {
        return String.format("type(%s) %s:%d", getType(), getChromosome(), getPosition());
    }
}
