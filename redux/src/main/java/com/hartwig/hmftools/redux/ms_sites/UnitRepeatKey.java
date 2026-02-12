package com.hartwig.hmftools.redux.ms_sites;

public class UnitRepeatKey
{
    public final UnitKey Key;
    public final int NumRepeats;

    public UnitRepeatKey(final UnitKey unitKey, final int numRepeats)
    {
        Key = unitKey;
        NumRepeats = numRepeats;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
        {
            return true;
        }
        if(!(o instanceof UnitRepeatKey))
        {
            return false;
        }

        final UnitRepeatKey that = (UnitRepeatKey) o;

        if(NumRepeats != that.NumRepeats)
        {
            return false;
        }
        return Key == that.Key;
    }

    @Override
    public int hashCode()
    {
        int result = Key.hashCode();
        result = 31 * result + NumRepeats;
        return result;
    }

    @Override
    public String toString()
    {
        return Key.getUnitKey() + " x " + NumRepeats;
    }
}
