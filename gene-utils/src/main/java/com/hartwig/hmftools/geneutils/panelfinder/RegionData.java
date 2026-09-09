package com.hartwig.hmftools.geneutils.panelfinder;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.mappability.RegionQuality;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.HighDepthRegion;

public class RegionData extends ChrBaseRegion
{
    private final List<HighDepthRegion> mHighDepths;
    private final List<GeneExonData> mGeneExons;
    private final List<PanelData> mPanelRegions;
    private final List<RegionQuality> mMappabilityScores;

    private String mPanelGeneName;
    private String mClosestGeneInfo;

    public RegionData(final ChrBaseRegion region)
    {
        super(region.Chromosome, region.start(), region.end());

        mHighDepths = Lists.newArrayList();
        mGeneExons = Lists.newArrayList();
        mPanelRegions = Lists.newArrayList();
        mMappabilityScores = Lists.newArrayList();
        mPanelGeneName = "";
        mClosestGeneInfo = "";
    }

    public List<HighDepthRegion> highDepths() { return mHighDepths; }
    public List<GeneExonData> geneExons() { return mGeneExons; }
    public List<PanelData> panelRegions() { return mPanelRegions; }

    public boolean panelRelated() { return !mPanelRegions.isEmpty() || !mPanelGeneName.isEmpty(); }

    public List<RegionQuality> mappabilityScores() { return mMappabilityScores; }

    public void setPanelGene(final String name) { mPanelGeneName = name; }
    public String panelGeneName() { return mPanelGeneName; }

    public void addHighDepth(final HighDepthRegion highDepth)
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

    public String closestGeneInfo() { return mClosestGeneInfo; }
    public void setClosestGeneInfo(final String info) { mClosestGeneInfo = info; }

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

        // limit probe qual regions to sections overlapping the region
        for(RegionQuality regionQuality : mMappabilityScores)
        {
            int regionOverlap = min(end(), regionQuality.end()) - max(start(), regionQuality.start()) + 1;
            totalQuality += regionQuality.Quality * regionOverlap;
            totalBases += regionOverlap;
        }

        return totalBases > 0 ? totalQuality / totalBases : 0;
    }

    public String label()
    {
        if(!panelRelated())
            return "NEW";

        if(mPanelRegions.size() == 1)
            return "EXISTING_" + mPanelRegions.get(0).Label;

        if(mPanelRegions.isEmpty() && !mPanelGeneName.isEmpty())
            return "NEW_" + mPanelGeneName;

        return format("EXISTING_%s_MULTI", mPanelRegions.get(0).Label);
    }

    public String toString()
    {
        return format("%s: highDepth(%s) genes(%s) panels(%s)",
                super.toString(), toString(mHighDepths), GeneExonData.toString(mGeneExons), PanelData.toString(mPanelRegions));
    }

    public static String toString(final List<HighDepthRegion> highDepths)
    {
        return highDepths.stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM));
    }

}
