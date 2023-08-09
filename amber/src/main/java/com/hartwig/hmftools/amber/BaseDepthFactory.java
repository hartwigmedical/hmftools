package com.hartwig.hmftools.amber;

import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.samtools.SamRecordUtils;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class BaseDepthFactory
{
    private final int mMinBaseQuality;

    public BaseDepthFactory(final int minBaseQuality)
    {
        mMinBaseQuality = minBaseQuality;
    }

    public void addEvidence(final BaseDepth evidence, final SAMRecord samRecord)
    {
        int quality = getBaseQuality(evidence, samRecord);
        if(quality >= mMinBaseQuality)
        {
            ++evidence.ReadDepth;

            int bafPosition = evidence.position();
            int readPosition = samRecord.getReadPositionAtReferencePosition(bafPosition);
            if(readPosition != 0)
            {
                if(!isIndel(bafPosition, readPosition, samRecord))
                {
                    char baseChar = samRecord.getReadString().charAt(readPosition - 1);

                    if(evidence.equalsRef(baseChar))
                    {
                        ++evidence.RefSupport;
                    }
                    else if(evidence.equalsAlt(baseChar))
                    {
                        ++evidence.AltSupport;
                    }
                }
                else
                {
                    ++evidence.IndelCount;
                }
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

    public static int getBaseQuality(@NotNull final GenomePosition position, @NotNull final SAMRecord samRecord)
    {
        // Get quality of base after del if necessary
        for(int pos = position.position(); pos <= samRecord.getAlignmentEnd(); pos++)
        {
            int readPosition = samRecord.getReadPositionAtReferencePosition(pos);
            if(readPosition != 0)
            {
                return SamRecordUtils.getBaseQuality(samRecord, readPosition);
            }
        }

        return 0;
    }

    public static BaseDepth copyBaseDepth(final BaseDepth pos)
    {
        return new BaseDepth(pos.Chromosome, pos.Position, pos.ref(), pos.alt());
    }

    public static BaseDepth fromAmberSite(final AmberSite site)
    {
        return new BaseDepth(site.chromosome(), site.position(), site.ref(), site.alt());
    }
}
