package com.hartwig.hmftools.lilac.seq;

import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class HlaSequence
{
    public final HlaAllele Allele;
    private String mRawSequence;

    public final HlaSequence copyWithAdditionalSequence(final String additionalSequence)
    {
        return new HlaSequence(Allele, mRawSequence + additionalSequence);
    }

    public void appendSequence(final String additionalSequence)
    {
        mRawSequence += additionalSequence;
    }

    public final String getRawSequence() { return mRawSequence; }

    public HlaSequence(final HlaAllele allele, final String rawSequence)
    {
        Allele = allele;
        mRawSequence = rawSequence;
    }

    public String toString()
    {
        return "HlaSequence(allele=" + Allele + ", rawSequence=" + mRawSequence + ")";
    }
}
