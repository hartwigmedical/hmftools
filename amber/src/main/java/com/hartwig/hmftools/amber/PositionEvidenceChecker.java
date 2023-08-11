package com.hartwig.hmftools.amber;

import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.samtools.SamRecordUtils;

import htsjdk.samtools.SAMRecord;

public class PositionEvidenceChecker
{
    private final int mMinBaseQuality;

    public PositionEvidenceChecker(final int minBaseQuality)
    {
        mMinBaseQuality = minBaseQuality;
    }

    public void addEvidence(final PositionEvidence posEvidence, final SAMRecord samRecord)
    {
        int baseQuality = getBaseQuality(posEvidence.Position, samRecord);

        if(baseQuality < mMinBaseQuality)
            return;

        ++posEvidence.ReadDepth;

        int bafPosition = posEvidence.position();
        int readPosition = samRecord.getReadPositionAtReferencePosition(bafPosition);
        if(readPosition != 0)
        {
            if(!isIndel(bafPosition, readPosition, samRecord))
            {
                char baseChar = samRecord.getReadString().charAt(readPosition - 1);

                if(posEvidence.equalsRef(baseChar))
                {
                    ++posEvidence.RefSupport;
                }
                else if(posEvidence.equalsAlt(baseChar))
                {
                    ++posEvidence.AltSupport;
                    posEvidence.AltQuality += baseQuality;
                }
            }
            else
            {
                ++posEvidence.IndelCount;
            }
        }
    }

    public static boolean isIndel(int bafPosition, int readPosition, final SAMRecord samRecord)
    {
        if(samRecord.getAlignmentEnd() > bafPosition)
        {
            // Delete?
            if(samRecord.getReadPositionAtReferencePosition(bafPosition + 1) == 0)
            {
                return true;
            }

            // Insert?
            return samRecord.getReferencePositionAtReadPosition(readPosition + 1) != bafPosition + 1;
        }

        return false;
    }

    public static int getBaseQuality(final int position, final SAMRecord samRecord)
    {
        // Get quality of base after del if necessary
        for(int pos = position; pos <= samRecord.getAlignmentEnd(); pos++)
        {
            int readPosition = samRecord.getReadPositionAtReferencePosition(pos);
            if(readPosition != 0)
            {
                return SamRecordUtils.getBaseQuality(samRecord, readPosition);
            }
        }

        return 0;
    }

    public static PositionEvidence fromAmberSite(final AmberSite site)
    {
        return new PositionEvidence(site.chromosome(), site.position(), site.ref(), site.alt());
    }
}
