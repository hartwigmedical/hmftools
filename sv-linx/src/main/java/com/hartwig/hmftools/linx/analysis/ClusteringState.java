package com.hartwig.hmftools.linx.analysis;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.linx.cn.HomLossEvent;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.SvBreakend;

public class ClusteringState
{
    private int mNextClusterId;

    private final Map<String, List<SvBreakend>> mChrBreakendMap; // every breakend on a chromosome, ordered by ascending position
    private List<LohEvent> mLohEventList;
    private List<HomLossEvent> mHomLossList;

    private long mDelCutoffLength;
    private long mDupCutoffLength;

    public static String CLUSTER_REASON_PROXIMITY = "Prox";
    public static String CLUSTER_REASON_LOH = "LOH";
    public static String CLUSTER_REASON_HOM_LOSS = "HomLoss";
    public static String CLUSTER_REASON_COMMON_ARMS = "ComArm";
    public static String CLUSTER_REASON_FOLDBACKS = "Foldback";
    public static String CLUSTER_REASON_SOLO_SINGLE = "Single";
    public static String CLUSTER_REASON_LOOSE_OVERLAP = "LooseOverlap";
    public static String CLUSTER_REASON_LOH_CHAIN = "LOHChain";
    public static String CLUSTER_REASON_NET_ARM_END_PLOIDY = "ArmEndPloidy";
    public static String CLUSTER_REASON_BE_PLOIDY_DROP = "BEPloidy";
    public static String CLUSTER_REASON_LONG_DEL_DUP_OR_INV = "DelDupInv";
    public static String CLUSTER_REASON_SATELLITE_SGL = "SatelliteSgl";

    public ClusteringState()
    {
        mChrBreakendMap = new HashMap();
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

    public long getDupCutoffLength() { return mDelCutoffLength; }
    public long getDelCutoffLength() { return mDupCutoffLength; }

    public void reset()
    {
        mNextClusterId = 0;
        mDelCutoffLength = 0;
        mDupCutoffLength = 0;
        mChrBreakendMap.clear();
    }

    public void setCutoffLengths(long delLength, long dupLength)
    {
        mDelCutoffLength = delLength;
        mDupCutoffLength = dupLength;
    }

    public void clearCopyNumberData()
    {
        if(mLohEventList != null)
        {
            mLohEventList.forEach(x -> x.clearBreakends());
        }

        if(mHomLossList != null)
        {
            mHomLossList.forEach(x -> x.clearBreakends());
        }
    }
}
