package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public class DiscordantGroup
{
    private final String mChromosome;
    private final Set<Integer> mJunctionPositions; // for info sake only
    private final Orientation mOrientation;
    private final List<Read> mReads;

    private final ChrBaseRegion mRemoteRegion;
    private final Orientation mRemoteOrientation;

    private int mMinAlignedPosition;
    private int mMaxAlignedPosition;

    private PhaseGroup mPhaseGroup;

    public DiscordantGroup(final Junction junction, final Read read)
    {
        mJunctionPositions = Sets.newHashSet(junction.Position);
        mChromosome = junction.Chromosome;
        mOrientation = junction.Orient;

        mReads = Lists.newArrayList(read);
        mRemoteRegion = new ChrBaseRegion(read.mateChromosome(), read.mateAlignmentStart(), read.mateAlignmentEnd());
        mRemoteOrientation = read.mateOrientation();

        mMinAlignedPosition = read.alignmentStart();
        mMaxAlignedPosition = read.alignmentEnd();
        mPhaseGroup = null;
    }

    public List<Read> reads() { return mReads; }

    public String chromosome() { return mChromosome; }
    public Orientation orientation() { return mOrientation; }

    public int minAlignedPosition() { return mMinAlignedPosition; }
    public int maxAlignedPosition() { return mMaxAlignedPosition; }

    public ChrBaseRegion remoteRegion() { return mRemoteRegion; }
    public Orientation remoteOrientation() { return mRemoteOrientation; }

    public PhaseGroup phaseGroup() { return mPhaseGroup; }
    public void setPhaseGroup(final PhaseGroup phaseGroup) { mPhaseGroup = phaseGroup; }

    public boolean matches(final Read read)
    {
        if(read.orientation() != mOrientation)
            return false;

        if(read.mateOrientation().isForward() != (mRemoteOrientation.isForward()))
            return false;

        // no buffer applied since will rely on a final merge for this
        if(!positionsOverlap(mMinAlignedPosition, mMaxAlignedPosition, read.alignmentStart(), read.alignmentEnd()))
            return false;

        return mRemoteRegion.overlaps(read.mateChromosome(), read.mateAlignmentStart(), read.mateAlignmentEnd());
    }

    public void addRead(final Read read)
    {
        mReads.add(read);

        mMinAlignedPosition = min(mMinAlignedPosition, read.alignmentStart());
        mMaxAlignedPosition = max(mMaxAlignedPosition, read.alignmentEnd());

        mRemoteRegion.setStart(min(mRemoteRegion.start(), read.mateAlignmentStart()));
        mRemoteRegion.setEnd(max(mRemoteRegion.end(), read.mateAlignmentEnd()));
    }

    public boolean hasRead(final Read read) { return mReads.stream().anyMatch(x -> x == read); }
    public boolean hasFragment(final String readId)
    {
        return mReads.stream().anyMatch(x -> x.id().equals(readId));
    }

    public boolean matches(final DiscordantGroup other)
    {
        if(other.orientation() != mOrientation)
            return false;

        if(other.remoteOrientation() != mRemoteOrientation)
            return false;


        if(!positionsOverlap(mMinAlignedPosition, mMaxAlignedPosition, other.minAlignedPosition(), other.maxAlignedPosition()))
            return false;

        return mRemoteRegion.overlaps(remoteRegion());
    }

    public String toString()
    {
        return format("range(%s:%d-%d) orient(%d) reads(%d) remote(%s orient=%d)",
                mChromosome, mMinAlignedPosition, mMaxAlignedPosition, mOrientation, mReads.size(), mRemoteRegion, mRemoteOrientation);
    }
}
