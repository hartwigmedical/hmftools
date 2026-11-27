package com.hartwig.hmftools.tars.liftback;

import com.hartwig.hmftools.common.codon.Nucleotides;

// Built per record on one worker thread, so the lazily filled reverse complement needs no synchronisation.
public final class ReadBases
{
    private final byte[] mBases;
    private byte[] mReverseComplement;

    public ReadBases(final byte[] bases)
    {
        mBases = bases;
    }

    public byte[] getBases() { return mBases; }

    public byte[] getReverseComplement()
    {
        if(mReverseComplement == null)
        {
            mReverseComplement = Nucleotides.reverseComplementBases(mBases);
        }
        return mReverseComplement;
    }
}
