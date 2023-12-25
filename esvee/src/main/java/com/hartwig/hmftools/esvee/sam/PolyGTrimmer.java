package com.hartwig.hmftools.esvee.sam;

import com.hartwig.hmftools.esvee.models.MutableRecord;
import com.hartwig.hmftools.esvee.util.Counter;

public class PolyGTrimmer
{
    private final Counter mRecordsTrimmed = new Counter("PolyG Trimmed");

    // This many G/Cs, or more
    private final int mPolyGCountThreshold;

    public PolyGTrimmer(final int polyGCountThreshold)
    {
        mPolyGCountThreshold = polyGCountThreshold;
    }

    public MutableRecord trimPolyG(final MutableRecord record)
    {
        if (record.isUnmapped() || !record.isPairedRead() || record.getLength() < mPolyGCountThreshold + 1)
            return record;

        if (record.isPositiveStrand())
        {
            int trailingGCount = 0;
            for(int i = record.getLength() - 1; i >= 0; i--)
                if(record.getBases()[i] == 'G')
                    trailingGCount++;
                else
                    break;
            if (trailingGCount >= mPolyGCountThreshold)
            {
                mRecordsTrimmed.add(1);
                return record.trimRight(trailingGCount);
            }
        }
        else
        {
            int leadingCCount = 0;
            for(int i = 0; i < record.getLength(); i++)
                if(record.getBases()[i] == 'C')
                    leadingCCount++;
                else
                    break;
            if (leadingCCount >= mPolyGCountThreshold)
            {
                mRecordsTrimmed.add(1);
                return record.trimLeft(leadingCCount);
            }
        }

        return record;
    }
}
