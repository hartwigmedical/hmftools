package com.hartwig.hmftools.linx.neoepitope;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RefGenomeSource implements RefGenomeInterface
{
    private final IndexedFastaSequenceFile mRefGenome;

    public RefGenomeSource(final IndexedFastaSequenceFile refGenome)
    {
        mRefGenome = refGenome;
    }

    @Override
    public String getBaseString(final String chromosome, int posStart, int posEnd)
    {
        return mRefGenome.getSubsequenceAt(chromosome, posStart, posEnd).getBaseString();
    }

}
