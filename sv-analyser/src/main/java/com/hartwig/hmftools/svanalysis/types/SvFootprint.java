package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.max;
import static java.lang.Math.min;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;

import java.util.List;

public class SvFootprint {

    private int mFootprintId;

    private List<SvClusterData> mSVs;
    private List<SvClusterData> mSpanningSVs;
    private List<SvFootprint> mLinkedFootprints;

    private final String mChromosome;
    private final String mArm;
    private long mStartPos;
    private long mEndPos;


    public SvFootprint(final int footprintId, final String chr, final String arm)
    {
        mFootprintId = footprintId;

        mChromosome = chr;
        mArm = arm;
        mStartPos = 0;
        mEndPos = 0;
        mSVs = Lists.newArrayList();
        mSpanningSVs = Lists.newArrayList();
        mLinkedFootprints = Lists.newArrayList();
    }

    public int getId() { return mFootprintId; }

    public final String posId() {
        return String.format("%d: %s_%s %d:%d", mFootprintId, mChromosome, mArm, mStartPos, mEndPos);
    }

    public final String chromosome() { return mChromosome; }
    public final String arm() { return mArm; }
    public long posStart() { return mStartPos; }
    public long posEnd() { return mEndPos; }

    public List<SvClusterData> getSVs() { return mSVs; }
    public List<SvClusterData> getSpanningSVs() { return mSpanningSVs; }
    public int getCount(boolean includeSpans) { return includeSpans ? mSVs.size() + mSpanningSVs.size() : mSVs.size(); }

    public void addVariant(final SvClusterData var, boolean isSpanning)
    {
        if(isSpanning) {
            mSpanningSVs.add(var);
        }
        else {
            mSVs.add(var);

            // widen the footprint range if required
            mStartPos = mStartPos == 0 ? var.position(true) : min(mStartPos, var.position(true));
            mEndPos = max(mEndPos, var.position(false));
        }
    }

    public void addLinkedFootprint(SvFootprint other) { mLinkedFootprints.add(other); }
    public final List<SvFootprint> getLinkedFootprints() { return mLinkedFootprints; }
}
