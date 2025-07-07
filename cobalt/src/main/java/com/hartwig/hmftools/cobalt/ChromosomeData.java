package com.hartwig.hmftools.cobalt;

public class ChromosomeData
{
    public final String Name;
    public final int Length;

    public ChromosomeData(String name, int length)
    {
        Name = name;
        Length = length;
    }

    @Override
    public String toString() { return Name; }
}
