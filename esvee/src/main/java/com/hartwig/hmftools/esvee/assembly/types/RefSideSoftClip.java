package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MAX_JUNC_POS_DIFF;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PHASED_ASSEMBLY_MIN_TI;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PROXIMATE_REF_SIDE_SOFT_CLIPS;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.REF_SIDE_MIN_SOFT_CLIP_LENGTH;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public class RefSideSoftClip
{
    public final int Position;
    public final Orientation Orient;

    private final List<String> mReadIds;
    private int mMaxLength;
    private boolean mMatchedOriginal;

    public RefSideSoftClip(final int position, final Orientation orientation, final Read read, final int readSoftClipLength)
    {
        Position = position;
        Orient = orientation;
        mMaxLength = 0;
        mReadIds = Lists.newArrayList();
        mMatchedOriginal = false;

        addRead(read, readSoftClipLength);
    }

    public List<String> readIds() { return mReadIds; }
    public int maxLength() { return mMaxLength; }
    public int readCount() { return mReadIds.size(); }

    public void markMatchesOriginal() { mMatchedOriginal = true; }
    public boolean matchesOriginal() { return mMatchedOriginal; }

    public void addRead(final Read read, final int readSoftClipLength)
    {
        mReadIds.add(read.id());
        mMaxLength = max(mMaxLength, readSoftClipLength);
    }

    public String toString() { return format("%d len(%d) reads(%d) %s",
            Position, mMaxLength, mReadIds.size(), mMatchedOriginal ? "matchOrig" : ""); }

    public static boolean checkAddRefSideSoftClip(final List<RefSideSoftClip> refSideSoftClips, final Junction junction, final Read read)
    {
        // collect reads which have soft-clips suggesting a short TI section ie a soft-clip on the other side the min TI length from the junction
        int refSideSoftClipPosition;
        int refSideSoftClipLength;

        if(junction.isForward())
        {
            refSideSoftClipPosition = read.alignmentStart();

            if(refSideSoftClipPosition > junction.Position - PHASED_ASSEMBLY_MIN_TI)
                return false;

            refSideSoftClipLength = read.leftClipLength();
        }
        else
        {
            refSideSoftClipPosition = read.alignmentEnd();

            if(refSideSoftClipPosition < junction.Position + PHASED_ASSEMBLY_MIN_TI)
                return false;

            refSideSoftClipLength = read.rightClipLength();
        }

        if(refSideSoftClipLength < ASSEMBLY_MAX_JUNC_POS_DIFF)
            return false;

        int softClipPosition = refSideSoftClipPosition;
        RefSideSoftClip existing = refSideSoftClips.stream().filter(x -> x.Position == softClipPosition).findFirst().orElse(null);

        if(existing == null)
        {
            refSideSoftClips.add(new RefSideSoftClip(softClipPosition, junction.Orient.opposite(), read, refSideSoftClipLength));
        }
        else
        {
            existing.addRead(read, refSideSoftClipLength);
        }

        return true;
    }

    public boolean hasProximateMatch(int otherRefPosition)
    {
        return abs(Position - otherRefPosition) <= PROXIMATE_REF_SIDE_SOFT_CLIPS;
    }

    public static void purgeRefSideSoftClips(final List<RefSideSoftClip> refSideSoftClips, int nonSoftClipRefPosition)
    {
        purgeRefSideSoftClips(refSideSoftClips, ASSEMBLY_MIN_READ_SUPPORT, REF_SIDE_MIN_SOFT_CLIP_LENGTH, nonSoftClipRefPosition);
    }

    public static void purgeRefSideSoftClips(
            final List<RefSideSoftClip> refSideSoftClips, int minCount, int minLength, int nonSoftClipRefPosition)
    {
        if(refSideSoftClips.isEmpty())
            return;

        // drop if not sufficient support or matches the original assembly's ref extension position anyway
        // or is close to it - where say homology causes aligned bases to match the soft-clip, or is past the non-soft-clipped position
        int index = 0;
        RefSideSoftClip matching = null;

        while(index < refSideSoftClips.size())
        {
            RefSideSoftClip refSideSoftClip = refSideSoftClips.get(index);

            boolean isPositionMatchOrClose = refSideSoftClip.hasProximateMatch(nonSoftClipRefPosition);
            boolean insufficientSupportOrLength = refSideSoftClip.readCount() < minCount || refSideSoftClip.maxLength() < minLength;
            boolean invalidPosition = refSideSoftClip.Orient.isForward() ?
                    refSideSoftClip.Position > nonSoftClipRefPosition : refSideSoftClip.Position < nonSoftClipRefPosition;

            if(insufficientSupportOrLength || isPositionMatchOrClose || invalidPosition)
            {
                refSideSoftClips.remove(index);

                if(isPositionMatchOrClose)
                {
                    if(refSideSoftClip.Position == nonSoftClipRefPosition)
                        matching = refSideSoftClip;
                    else if(matching == null || matching.readCount() < refSideSoftClip.readCount())
                        matching = refSideSoftClip;
                }
            }
            else
            {
                ++index;
            }
        }

        // retain info about any matching soft-clip
        if(matching != null)
        {
            matching.markMatchesOriginal();
            refSideSoftClips.add(matching);
        }
    }
}
