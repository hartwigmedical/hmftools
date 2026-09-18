package com.hartwig.hmftools.virusdetect;

// A group of viruses at the level of taxonomy granularity which matters to us. Detection resolves one representative
// contig per group, so the group is the unit the output reports on.
public record OncologyGroup(String name)
{
    public OncologyGroup
    {
        if(name.isEmpty())
        {
            throw new IllegalArgumentException("oncology group name is empty");
        }
    }

    @Override
    public String toString()
    {
        return name;
    }
}
