package com.hartwig.hmftools.tars.liftback.features;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.tars.liftback.LiftedAlignment;

import htsjdk.samtools.SAMRecord;

// Built per record on one worker thread, so the lazily filled reverse complement needs no synchronisation.
final class ReadBases
{
    private final byte[] mForward;
    private final boolean mRecordForward;
    private byte[] mReverse;

    private ReadBases(final byte[] forward, final boolean recordForward)
    {
        mForward = forward;
        mRecordForward = recordForward;
    }

    // null when the record carries no sequence, which leaves nothing to score or gate against
    static ReadBases of(final SAMRecord record)
    {
        byte[] forward = record.getReadBases();
        if(forward == null || forward.length == 0)
        {
            return null;
        }

        return new ReadBases(forward, !record.getReadNegativeStrandFlag());
    }

    byte[] forAlignment(final LiftedAlignment alignment)
    {
        if(alignment.ForwardStrand == mRecordForward)
        {
            return mForward;
        }

        if(mReverse == null)
        {
            mReverse = Nucleotides.reverseComplementBases(mForward);
        }
        return mReverse;
    }
}
