package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.common.amber.BaseDepthFactory.getBaseQuality;
import static com.hartwig.hmftools.common.amber.BaseDepthFactory.indel;

import com.hartwig.hmftools.common.amber.BaseDepth;

import htsjdk.samtools.SAMRecord;

class TumorBAFFactory
{
    private final int mMinBaseQuality;

    TumorBAFFactory(final int minBaseQuality)
    {
        mMinBaseQuality = minBaseQuality;
    }

    public static TumorBAF create(final BaseDepth normal)
    {
        TumorBAF tumorBAF = new TumorBAF(normal.chromosome(), normal.position(), normal.ref().toString(), normal.alt().toString());
        tumorBAF.NormalReadDepth = normal.readDepth();
        tumorBAF.NormalRefSupport = normal.refSupport();
        tumorBAF.NormalAltSupport = normal.altSupport();
        return tumorBAF;
    }

    void addEvidence(final TumorBAF evidence, final SAMRecord samRecord)
    {
        int quality = getBaseQuality(evidence, samRecord);
        if(quality >= mMinBaseQuality)
        {
            evidence.TumorReadDepth = evidence.TumorReadDepth + 1;

            int bafPosition = (int) evidence.position();
            int readPosition = samRecord.getReadPositionAtReferencePosition(bafPosition);
            if(readPosition != 0)
            {
                if(!indel(bafPosition, readPosition, samRecord))
                {
                    final String base = String.valueOf(samRecord.getReadString().charAt(readPosition - 1));
                    if(base.equals(evidence.Ref))
                    {
                        ++evidence.TumorRefSupport;
                    }
                    else if(base.equals(evidence.Alt))
                    {
                        ++evidence.TumorAltSupport;
                        evidence.TumorAltQuality += quality;
                    }
                }
                else
                {
                    ++evidence.TumorIndelCount;
                }
            }
        }
    }
}
