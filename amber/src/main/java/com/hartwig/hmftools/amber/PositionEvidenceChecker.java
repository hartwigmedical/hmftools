package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConfig.isUltima;

import java.util.List;

import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.sequencing.UltimaBamUtils;

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

        int bafPosition = posEvidence.position();
        int readIndex = read.getReadPositionAtReferencePosition(bafPosition) - 1; // is 1-based

        if(isUltima())
        {
            List<Integer> lowQualIndices = UltimaBamUtils.extractLowQualIndices(read);

            if(lowQualIndices != null && lowQualIndices.contains(readIndex))
            {
                ++posEvidence.SeqTechFiltered;
                filtered = true;
            }
        }

        ++posEvidence.ReadDepth;

        if(filtered)
            return;

        if(readIndex >= 0)
        {
            if(!isIndel(bafPosition, readIndex + 1, read)) // to preserve existing read index behaviour
            {
                char baseChar = read.getReadString().charAt(readIndex);

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

    private static boolean isIndel(int bafPosition, int readIndex, final SAMRecord samRecord)
    {
        if(samRecord.getAlignmentEnd() > bafPosition)
        {
            // Delete?
            if(samRecord.getReadPositionAtReferencePosition(bafPosition + 1) == 0)
            {
                return true;
            }

            // Insert?
            return samRecord.getReferencePositionAtReadPosition(readIndex + 1) != bafPosition + 1;
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
