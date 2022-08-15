package com.hartwig.hmftools.linx.drivers;

import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.linx.cn.SvCNData;

public class GeneCopyNumberRegion
{
    public final String GeneName;
    public final String TransName;
    public final int RegionStart;
    public final int RegionEnd;
    public final double MinCopyNumber;
    public final boolean IsGermline;

    public GeneCopyNumberRegion(final GeneCopyNumber geneCN)
    {
        GeneName = geneCN.geneName();
        TransName = geneCN.transName();
        RegionStart = (int)geneCN.minRegionStart();
        RegionEnd = (int)geneCN.minRegionEnd();
        MinCopyNumber = geneCN.minCopyNumber();
        IsGermline = false; // TODO: check new logic geneCN.germlineHet2HomRegions() > 0 || geneCN.germlineHomRegions() > 0;
    }

    public GeneCopyNumberRegion(final String geneName, final String transName, int regionStart, int regionEnd, double minCopyNumber)
    {
        GeneName = geneName;
        TransName = transName;
        RegionStart = regionStart;
        RegionEnd = regionEnd;
        MinCopyNumber = minCopyNumber;
        IsGermline = false;
    }

    public static GeneCopyNumberRegion calcGeneCopyNumberRegion(final TranscriptData transData, final List<SvCNData> copyNumberData)
    {
        SvCNData minRegion = null;

        for(final ExonData exon : transData.exons())
        {
            for(SvCNData cnData : copyNumberData)
            {
                if(cnData.EndPos < exon.Start)
                    continue;

                if(cnData.StartPos > exon.End)
                    break;

                if(minRegion == null || cnData.CopyNumber < minRegion.CopyNumber)
                {
                    minRegion = cnData;
                }
            }
        }

        return minRegion != null ? new GeneCopyNumberRegion(
                transData.GeneId, transData.TransName, minRegion.StartPos, minRegion.EndPos, minRegion.CopyNumber) : null;
    }

}
