package com.hartwig.hmftools.geneutils.panelfinder;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class RegionData extends ChrBaseRegion
{
    private final List<HighDepthData> mHighDepths;
    private final List<GeneExonData> mGeneExons;
    private final List<PanelData> mPanelRegions;

    private final List<Double> mMappabilityScores;

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
    public List<Double> mappabilityScores() { return mMappabilityScores; }

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

    public String toString()
    {
        return String.format("%s: highDepth(%s) genes(%s) panels(%s)",
                super.toString(), HighDepthData.toString(mHighDepths), GeneExonData.toString(mGeneExons), PanelData.toString(mPanelRegions));
    }

    /*
    public static void mergeRegion(final List<RegionData> regions, final RegionData newRegion)
    {
        // insert regions in ascending order by position
        // merge any overlapping regions
        int index = 0;

        while(index < regions.size())
        {
            RegionData region = regions.get(index);

            if(newRegion.start() > region.end())
            {
                ++index;
                continue;
            }

            if(region.start() > newRegion.end())
                break;

            if(newRegion.matches(region))
                return;

            // handle merges
            int startPosition = min(region.start(), newRegion.start());
            region.setStart(startPosition);

            int endPosition = max(region.end(), newRegion.end());
            region.setEnd(endPosition);

            ++index;

            while(index < regions.size())
            {
                RegionData nextRegion = regions.get(index);

                if(nextRegion.start() > region.end())
                    break;

                endPosition = max(region.end(), nextRegion.end());
                region.setEnd(endPosition);
                regions.remove(index);
            }

            return;
        }

        regions.add(index, newRegion);
    }

    public static void integrateRegion(final List<RegionData> regions, final RegionData newRegion)
    {
        // split these new regions if they overlap an existing coding region
        int index = 0;

        while(index < regions.size())
        {
            RegionData region = regions.get(index);

            if(newRegion.start() > region.end())
            {
                ++index;
                continue;
            }

            if(region.start() > newRegion.end())
                break;

            if(newRegion.matches(region))
                return;

            if(newRegion.start() < region.start())
            {
                RegionData preRegion = new RegionData(
                        newRegion.GeneName,
                        new ChrBaseRegion(newRegion.Chromosome, newRegion.start(), region.start() - 1),
                        newRegion.ExonRank, newRegion.Type);

                preRegion.setExtraInfo(newRegion.getExtraInfo());

                regions.add(index, preRegion);
                ++index; // for the additional insert

                // adjust for remaining segment
                if(newRegion.end() <= region.end())
                    return;

                newRegion.setStart(region.end() + 1);
            }
            else
            {
                newRegion.setStart(region.end() + 1);
            }

            ++index;
        }

        regions.add(index, newRegion);
    }

    public static boolean validate(final List<RegionData> regions)
    {
        for(int i = 0; i < regions.size() - 1; ++i)
        {
            RegionData region = regions.get(i);
            RegionData nextRegion = regions.get(i + 1);

            if(nextRegion.start() <= region.end())
                return false;
        }

        return true;
    }
    */
}
