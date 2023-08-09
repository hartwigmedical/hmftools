package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.BaseDepthFactory.getBaseQuality;
import static com.hartwig.hmftools.amber.BaseDepthFactory.isIndel;

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
        TumorBAF tumorBAF = new TumorBAF(normal.chromosome(), normal.position(), normal.ref(), normal.alt());
        tumorBAF.NormalReadDepth = normal.ReadDepth;
        tumorBAF.NormalRefSupport = normal.RefSupport;
        tumorBAF.NormalAltSupport = normal.AltSupport;
        return tumorBAF;
    }

    void addEvidence(final TumorBAF evidence, final SAMRecord samRecord)
    {
        int quality = getBaseQuality(evidence, samRecord);
        if(quality >= mMinBaseQuality)
        {
            evidence.TumorReadDepth = evidence.TumorReadDepth + 1;

            int bafPosition = evidence.position();
            int readPosition = samRecord.getReadPositionAtReferencePosition(bafPosition);
            if(readPosition != 0)
            {
                if(!isIndel(bafPosition, readPosition, samRecord))
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
