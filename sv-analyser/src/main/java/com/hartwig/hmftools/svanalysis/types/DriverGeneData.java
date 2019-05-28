package com.hartwig.hmftools.svanalysis.types;

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
    private boolean mFullyMatched;
    private boolean mMissedLohSVs;

    public DriverGeneData(final DriverCatalog driverGene, final HmfTranscriptRegion region, final GeneCopyNumber geneCN)
    {
        DriverGene = driverGene;
        Region = region;
        GeneCN = geneCN;
        mMissedLohSVs = false;
        mFullyMatched = false;

        mSvBreakends = Lists.newArrayList();
        mSvInfo = Lists.newArrayList();
    }

    public void setMissedLohSVs(boolean toggle) { mMissedLohSVs = toggle; }
    public boolean missedLohSVs() { return mMissedLohSVs; }

    public void setFullyMatched(boolean toggle) { mFullyMatched = toggle; }
    public boolean fullyMatched() { return mFullyMatched; }

    public final List<SvBreakend> getSvBreakends() { return mSvBreakends; }
    public final List<String> getSvInfoList() { return mSvInfo; }

    public void addSvBreakend(final SvBreakend breakend, final String info)
    {
        mSvBreakends.add(breakend);
        mSvInfo.add(info);
    }

    public void addMatchInfo(final String info)
    {
        mSvInfo.add(info);
    }

}
