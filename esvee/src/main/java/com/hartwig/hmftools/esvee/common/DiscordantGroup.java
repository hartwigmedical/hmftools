package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.common.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.common.SupportType.INDEL;
import static com.hartwig.hmftools.esvee.common.SupportType.JUNCTION;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.read.ReadUtils;

public class DiscordantGroup
{
    private final Junction mJunction;
    private final List<Read> mReads;
    private final ChrBaseRegion mRemoteRegion;

    private int mMinAlignedPosition;
    private int mMaxAlignedPosition;

    private byte mBases[];
    private byte mBaseQuals[];

    private PhaseGroup mPhaseGroup;

    public DiscordantGroup(final Junction junction, final Read read)
    {
        mJunction = junction;
        mReads = Lists.newArrayList(read);
        mRemoteRegion = new ChrBaseRegion(read.mateChromosome(), read.mateAlignmentStart(), read.mateAlignmentEnd());

        mMinAlignedPosition = 0;
        mMaxAlignedPosition = 0;
        mBases = null;
        mBaseQuals = null;

        mPhaseGroup = null;
    }

    public Junction junction() { return mJunction; }
    public List<Read> reads() { return mReads; }
    public ChrBaseRegion remoteRegion() { return mRemoteRegion; }

    public int minAlignedPosition() { return mMinAlignedPosition; }
    public int maxAlignedPosition() { return mMaxAlignedPosition; }

    public PhaseGroup phaseGroup() { return mPhaseGroup; }
    public void setPhaseGroup(final PhaseGroup phaseGroup) { mPhaseGroup = phaseGroup; }

    public boolean matches(final Read read)
    {
        return mRemoteRegion.overlaps(read.mateChromosome(), read.mateAlignmentStart(), read.mateAlignmentEnd());
    }

    public void addRead(final Read read)
    {
        mReads.add(read);

        mRemoteRegion.setStart(min(mRemoteRegion.start(), read.mateAlignmentStart()));
        mRemoteRegion.setEnd(max(mRemoteRegion.end(), read.mateAlignmentEnd()));
    }

    public boolean hasRead(final Read read)
    {
        return read != null && mReads.stream().anyMatch(x -> x.matchesFragment(read));
    }

    public void buildAssembly()
    {
        // establish the boundaries from the reads but for now don't bother building a sequence
        // if were to, would be very similar to how RefBaseAssembly is done
        int minAlignedPosition = mJunction.Position;
        int maxAlignedPosition = mJunction.Position;

        for(Read read : mReads)
        {
            if(mJunction.isForward())
            {
                maxAlignedPosition = max(maxAlignedPosition, read.unclippedEnd());
            }
            else
            {
                minAlignedPosition = min(minAlignedPosition, read.unclippedStart());
            }
        }

        mMinAlignedPosition = minAlignedPosition;
        mMaxAlignedPosition = maxAlignedPosition;


        /*
        Read maxBaseQualRead = null;
        double maxAvgBaseQual = 0;

        int maxDistanceFromJunction = 0;

        // could sort and build from reads without gaps, but depends on how used in linking - likely not required
        for(Read read : mReads)
        {
            int readJunctionIndex = read.getReadIndexAtReferencePosition(mJunction.Position, true);

            int readJunctionDistance;

            if(mJunction.isForward())
            {
                maxAlignedPosition = max(maxAlignedPosition, read.unclippedEnd());
                readJunctionDistance = readJunctionIndex;
            }
            else
            {
                minAlignedPosition = min(minAlignedPosition, read.unclippedStart());
                readJunctionDistance = read.basesLength() - readJunctionIndex;
            }

            maxDistanceFromJunction = max(maxDistanceFromJunction, readJunctionDistance);

            double avgBaseQual = ReadUtils.avgBaseQuality(read);

            if(avgBaseQual > maxAvgBaseQual)
            {
                maxAvgBaseQual = avgBaseQual;
                maxBaseQualRead = read;
            }
        }

        for(Read read : mReads)
        {


        }

        // buildRepeatInfo();
        */
    }

    public String toString()
    {
        return format("reads(%d) remote(%s)", mReads.size(), mRemoteRegion);
    }
}
