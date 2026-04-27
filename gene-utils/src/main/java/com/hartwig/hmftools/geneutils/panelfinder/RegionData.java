package com.hartwig.hmftools.geneutils.panelfinder;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.mappability.RegionQuality;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class RegionData extends ChrBaseRegion
{
    private final List<HighDepthData> mHighDepths;
    private final List<GeneExonData> mGeneExons;
    private final List<PanelData> mPanelRegions;
    private final List<RegionQuality> mMappabilityScores;

    public RegionData(final ChrBaseRegion region)
    {
        super(region.Chromosome, region.start(), region.end());

        mHighDepths = Lists.newArrayList();
        mGeneExons = Lists.newArrayList();
        mPanelRegions = Lists.newArrayList();
        mMappabilityScores = Lists.newArrayList();
    }

    public List<HighDepthData> highDepths() { return mHighDepths; }
    public List<GeneExonData> geneExons() { return mGeneExons; }
    public List<PanelData> panelRegions() { return mPanelRegions; }
    public List<RegionQuality> mappabilityScores() { return mMappabilityScores; }

    public void addHighDepth(final HighDepthData highDepth)
    {
        mHighDepths.add(highDepth);

        setStart(min(start(), highDepth.start()));
        setEnd(max(end(), highDepth.end()));
    }

    public void addPanelData(final PanelData panelData)
    {
        mPanelRegions.add(panelData);

        setStart(min(start(), panelData.start()));
        setEnd(max(end(), panelData.end()));
    }

    public void addGeneExon(final GeneExonData geneExon)
    {
        mGeneExons.add(geneExon);
    }

    public void mergeRegion(final RegionData other)
    {
        setStart(min(start(), other.start()));
        setEnd(max(end(), other.end()));

        mHighDepths.addAll(other.highDepths());
        mGeneExons.addAll(other.geneExons());
        mPanelRegions.addAll(other.panelRegions());
    }

    public double meanMappability()
    {
        if(mMappabilityScores.isEmpty())
            return 0;

        int totalBases = 0;
        double totalQuality = 0;

        for(RegionQuality regionQuality : mMappabilityScores)
        {
            totalQuality += regionQuality.Quality * regionQuality.baseLength();
            totalBases += regionQuality.baseLength();
        }

        return totalBases > 0 ? totalQuality / totalBases : 0;
    }

    public String toString()
    {
        return String.format("%s: highDepth(%s) genes(%s) panels(%s)",
                super.toString(), HighDepthData.toString(mHighDepths), GeneExonData.toString(mGeneExons), PanelData.toString(mPanelRegions));
    }
}
