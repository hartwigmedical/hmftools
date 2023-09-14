package com.hartwig.hmftools.amber;

import static java.lang.Math.max;
import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;

public class RegionTask
{
    private final PositionEvidenceChecker mEvidenceChecker;
    public final ChrBaseRegion Region;

    private final List<PositionEvidence> mPositions;
    private int mCurrentIndex;
    private boolean mComplete;

    public RegionTask(final PositionEvidenceChecker evidenceChecker, final String chromosome, final PositionEvidence baseDepth)
    {
        mEvidenceChecker = evidenceChecker;
        Region = new ChrBaseRegion(chromosome, baseDepth.Position, baseDepth.Position);
        mPositions = Lists.newArrayList(baseDepth);
        mCurrentIndex = 0;
        mComplete = false;
    }

    public void addPosition(final PositionEvidence posEvidence)
    {
        mPositions.add(posEvidence);
        Region.setEnd(max(Region.end(), posEvidence.Position));
    }

    public void processRecord(final SAMRecord record)
    {
        int alignmentStart = record.getAlignmentStart();
        int alignmentEnd = record.getAlignmentEnd();

        int index = mCurrentIndex;
        for(; index < mPositions.size(); ++index)
        {
            PositionEvidence posEvidence = mPositions.get(index);

            if(alignmentStart > posEvidence.Position)
            {
                ++mCurrentIndex;
                continue;
            }

            if(alignmentEnd < posEvidence.Position)
                break;

            mEvidenceChecker.addEvidence(posEvidence, record);
        }

        if(mCurrentIndex >= mPositions.size())
            mComplete = true;
    }

    public boolean isComplete()
    {
        return mComplete;
    }

    public int positionCount()
    {
        return mPositions.size();
    }

    public String toString()
    {
        return format("region(%s) positions(%d) index(%d)", Region, mPositions.size(), mCurrentIndex);
    }
}
