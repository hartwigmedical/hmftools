package com.hartwig.hmftools.linx.analysis;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.cn.HomLossEvent;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.SvBreakend;

public class ClusteringState
{
    private int mNextClusterId;

    // every breakend on a chromosome, ordered by ascending position
    private final Map<String, List<SvBreakend>> mChrBreakendMap;

    private List<LohEvent> mLohEventList;
    private List<HomLossEvent> mHomLossList;

    private int mDelCutoffLength;
    private int mDupCutoffLength;

    public ClusteringState()
    {
        mChrBreakendMap = Maps.newHashMap();
        mLohEventList = null;
        mHomLossList = null;

        reset();
    }

    public final Map<String, List<SvBreakend>> getChrBreakendMap() { return mChrBreakendMap; }
    public final List<LohEvent> getLohEventList() { return mLohEventList; }
    public final List<HomLossEvent> getHomLossList() { return mHomLossList; }
    public int getNextClusterId() { return mNextClusterId++; }

    public void setSampleCnEventData(final List<LohEvent> lohEvents, List<HomLossEvent> homLossEvents)
    {
        mLohEventList = lohEvents;
        mHomLossList = homLossEvents;
    }

    public int getDupCutoffLength() { return mDelCutoffLength; }
    public int getDelCutoffLength() { return mDupCutoffLength; }

    public void reset()
    {
        mNextClusterId = 0;
        mDelCutoffLength = 0;
        mDupCutoffLength = 0;
        mChrBreakendMap.clear();
    }

    public void setCutoffLengths(int delLength, int dupLength)
    {
        mDelCutoffLength = delLength;
        mDupCutoffLength = dupLength;
    }
}
