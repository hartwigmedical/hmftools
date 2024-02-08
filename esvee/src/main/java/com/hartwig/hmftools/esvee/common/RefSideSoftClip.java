package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.SvConstants.PROXIMATE_REF_SIDE_SOFT_CLIPS;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.read.Read;

public class RefSideSoftClip
{
    public final int Position;
    public final byte Orientation;

    private final List<Read> mReads;
    private int mMaxLength;
    private boolean mMatchedOriginal;

    public RefSideSoftClip(final int position, final byte orientation, final Read read, final int readSoftClipLength)
    {
        Position = position;
        Orientation = orientation;
        mReads = Lists.newArrayList(read);
        mMaxLength = readSoftClipLength;
        mMatchedOriginal = false;
    }

    public List<Read> reads() { return mReads; }
    public int maxLength() { return mMaxLength; }
    public int readCount() { return mReads.size(); }

    public void markMatchesOriginal() { mMatchedOriginal = true; }
    public boolean matchesOriginal() { return mMatchedOriginal; }

    public void addRead(final Read read, final int readSoftClipLength)
    {
        mReads.add(read);
        mMaxLength = max(mMaxLength, readSoftClipLength);
    }

    public String toString() { return format("%d len(%d) reads(%d)", Position, mMaxLength, mReads.size()); }

    public static void purgeRefSideSoftClips(
            final List<RefSideSoftClip> refSideSoftClips, int minCount, int minLength, int nonSoftClipRefPosition)
    {
        if(refSideSoftClips.isEmpty())
            return;

        // drop if not sufficient support or matches the original assembly's ref extension position anyway
        // or is close to it - where say homology causes aligned bases to match the soft-clip
        int index = 0;
        RefSideSoftClip matching = null;

        while(index < refSideSoftClips.size())
        {
            RefSideSoftClip refSideSoftClip = refSideSoftClips.get(index);

            boolean isPositionMatchOrClose = abs(refSideSoftClip.Position - nonSoftClipRefPosition) <= PROXIMATE_REF_SIDE_SOFT_CLIPS;

            if(refSideSoftClip.readCount() < minCount || refSideSoftClip.maxLength() < minLength || isPositionMatchOrClose)
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
