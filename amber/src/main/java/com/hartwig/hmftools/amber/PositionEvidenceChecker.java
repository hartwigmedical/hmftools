package com.hartwig.hmftools.amber;

import com.hartwig.hmftools.common.amber.AmberSite;

import htsjdk.samtools.SAMRecord;

public class PositionEvidenceChecker
{
    private final int mMinMapQuality;
    private final int mMinBaseQuality;

    public PositionEvidenceChecker(final int minMapQuality, final int minBaseQuality)
    {
        mMinMapQuality = minMapQuality;
        mMinBaseQuality = minBaseQuality;
    }

    public void addEvidence(final PositionEvidence posEvidence, final SAMRecord read)
    {
        int baseQuality = getBaseQuality(posEvidence.Position, read);

        boolean filtered = false;

        if(read.getMappingQuality() < mMinMapQuality)
        {
            ++posEvidence.MapQualFiltered;
            filtered = true;
        }

        if(baseQuality < mMinBaseQuality)
        {
            ++posEvidence.BaseQualFiltered;
            filtered = true;
        }

        ++posEvidence.ReadDepth;

        if(filtered)
            return;

        int bafPosition = posEvidence.position();
        int readPosition = read.getReadPositionAtReferencePosition(bafPosition);
        if(readPosition != 0)
        {
            if(!isIndel(bafPosition, readPosition, read))
            {
                char baseChar = read.getReadString().charAt(readPosition - 1);

                if(posEvidence.equalsRef(baseChar))
                {
                    ++posEvidence.RefSupport;
                }
                else if(posEvidence.equalsAlt(baseChar))
                {
                    ++posEvidence.AltSupport;
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

    private static int getBaseQuality(final int position, final SAMRecord samRecord)
    {
        // get quality of base after del if necessary
        for(int pos = position; pos <= samRecord.getAlignmentEnd(); pos++)
        {
            int readPosition = samRecord.getReadPositionAtReferencePosition(pos);

            if(readPosition > 0)
            {
                return samRecord.getBaseQualities()[readPosition - 1];
            }
        }

        return 0;
    }

    public static PositionEvidence fromAmberSite(final AmberSite site)
    {
        return new PositionEvidence(site.chromosome(), site.position(), site.ref(), site.alt());
    }
}
