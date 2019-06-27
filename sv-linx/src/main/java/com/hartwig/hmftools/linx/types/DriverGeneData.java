package com.hartwig.hmftools.linx.types;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;

public class DriverGeneData
{
    public final DriverCatalog DriverGene;
    public final HmfTranscriptRegion Region;
    public final GeneCopyNumber GeneCN;

    private List<SvBreakend> mSvBreakends;
    private List<String> mSvInfo;
    private List<Integer> mLinkingIds;
    private boolean mFullyMatched;

    public DriverGeneData(final DriverCatalog driverGene, final HmfTranscriptRegion region, final GeneCopyNumber geneCN)
    {
        DriverGene = driverGene;
        Region = region;
        GeneCN = geneCN;
        mFullyMatched = false;

        mSvBreakends = Lists.newArrayList();
        mSvInfo = Lists.newArrayList();
        mLinkingIds = Lists.newArrayList();
    }

    public void setFullyMatched(boolean toggle) { mFullyMatched = toggle; }
    public boolean fullyMatched() { return mFullyMatched; }

    public final List<SvBreakend> getSvBreakends() { return mSvBreakends; }
    public final List<String> getSvInfoList() { return mSvInfo; }
    public final List<Integer> getLinkingIds() { return mLinkingIds; }

    public void addSvBreakend(final SvBreakend breakend, final String info)
    {
        mSvBreakends.add(breakend);
        mSvInfo.add(info);
    }

    public void addSvBreakendPair(final SvBreakend breakend1, final SvBreakend breakend2, final String info)
    {
        mSvBreakends.add(breakend1);
        mSvBreakends.add(breakend2);
        mSvInfo.add(info);
        mSvInfo.add(info);

        int linkingId = mLinkingIds.size() / 2;
        mLinkingIds.add(linkingId);
        mLinkingIds.add(linkingId);
    }

    public void addMatchInfo(final String info)
    {
        mSvInfo.add(info);
    }

}
