package com.hartwig.hmftools.common.genome.tiny;

import java.util.HashSet;
import java.util.Set;

public class TinyGenome extends SimpleTestGenome
{
    private final Set<ChromosomeSnippet> data = new HashSet<>();

    public void add(ChromosomeSnippet snippet)
    {
        data.add(snippet);
    }

    @Override
    public byte[] getBases(final String chromosome, final int posStart, final int posEnd)
    {
        ChromosomeSnippet relevantSnippet = data.stream()
                .filter(snippet -> snippet.containsData(chromosome, posStart, posEnd))
                .findFirst().orElseThrow();
        return relevantSnippet.getBases(posStart, posEnd);
    }

    @Override
    public String getBaseString(final String chromosome, final int posStart, final int posEnd)
    {
        return new String(getBases(chromosome, posStart, posEnd));
    }
}

