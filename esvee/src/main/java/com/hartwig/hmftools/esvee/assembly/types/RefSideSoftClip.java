package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MAX_JUNC_POS_DIFF;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_SPLIT_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PHASED_ASSEMBLY_MIN_TI;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PROXIMATE_REF_SIDE_SOFT_CLIPS;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.REF_SIDE_MIN_SOFT_CLIP_LENGTH;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public class RefSideSoftClip
{
    public final int Position;
    public final Orientation Orient;

    private final Set<String> mReadIds;
    private int mMaxLength;
    private boolean mMatchedOriginal;

    public RefSideSoftClip(final int position, final Orientation orientation, final Read read, final int readSoftClipLength)
    {
        Position = position;
        Orient = orientation;
        mMaxLength = 0;
        mReadIds = Sets.newHashSet();
        mMatchedOriginal = false;

        addRead(read, readSoftClipLength);
    }

    public Set<String> readIds() { return mReadIds; }
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
        // collect reads which have soft-clips suggesting a short TI section ie a soft-clip on the other side
        // of the min TI length from the junction
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

    public static void checkSupportVsRefSideSoftClip(final JunctionAssembly assembly)
    {
        // purge any junction support which extends beyond a consensus ref-side soft-clip
        List<RefSideSoftClip> refSideSoftClips = Lists.newArrayList();

        for(SupportRead read : assembly.support())
        {
            checkAddRefSideSoftClip(refSideSoftClips, assembly.junction(), read.cachedRead());
        }

        if(refSideSoftClips.isEmpty())
            return;

        Collections.sort(refSideSoftClips, Comparator.comparingInt(x -> -x.readCount()));

        RefSideSoftClip refSideSoftClip = refSideSoftClips.get(0);

        if(refSideSoftClip.readCount() < ASSEMBLY_MIN_READ_SUPPORT)
            return;

        int mainSoftClipCount = refSideSoftClip.readCount();
        int totalSoftClipCount = refSideSoftClips.stream().mapToInt(x -> x.readCount()).sum();
        int nonSoftClipCount = assembly.supportCount() - totalSoftClipCount;

        if(nonSoftClipCount > 0 && nonSoftClipCount < mainSoftClipCount && nonSoftClipCount < ASSEMBLY_SPLIT_MIN_READ_SUPPORT)
        {
            // as per the branching routine run during linking, require a minimum number of reads to keep both reads which soft-clip and
            // those which run past that point
            List<SupportRead> support = assembly.support();
            int index = 0;

            while(index < support.size())
            {
                SupportRead read = support.get(index);

                if(refSideSoftClips.stream().noneMatch(x -> x.readIds().contains(read.id())))
                {
                    support.remove(index);
                }
                else
                {
                    ++index;
                }
            }
        }
    }
}
