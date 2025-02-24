package com.hartwig.hmftools.pave.transval;

import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;

public interface TranscriptFilter
{
    boolean applies(TranscriptAminoAcids aminoAcids);
}
class Negation implements TranscriptFilter
{
    private final TranscriptFilter filter;

    public Negation(final TranscriptFilter filter)
    {
        this.filter = filter;
    }

    @Override
    public boolean applies(final TranscriptAminoAcids aminoAcids)
    {
        return !filter.applies(aminoAcids);
    }

    @Override
    public String toString()
    {
        return "Negation{" +
                "filter=" + filter +
                '}';
    }
}

class PassThroughFilter implements TranscriptFilter
{
    @Override
    public boolean applies(final TranscriptAminoAcids aminoAcids)
    {
        return true;
    }

    @Override
    public String toString()
    {
        return "PassThroughFilter{}";
    }
}
