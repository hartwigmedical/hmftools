package com.hartwig.hmftools.linx.fusion_likelihood;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.haveOverlap;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.PHASE_MAX;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.intAsType;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class PhaseRegionUtils
{
    private static final Logger LOGGER = LogManager.getLogger(PhaseRegionUtils.class);

    public static void checkAddCombinedGenePhaseRegion(final GenePhaseRegion regionToAdd, final List<GenePhaseRegion> regions)
    {
        // combine all regions into non-overlapping regions, allowing for mixed phasings

        // clone to avoid changing the contents of the original
        GenePhaseRegion tmp = GenePhaseRegion.from(regionToAdd);

        List<GenePhaseRegion> newRegions = Lists.newArrayList(tmp);

        // look for overlapping regions and combine or split them as required
        // split any new regions until they can be added without any further splits
        while(!newRegions.isEmpty())
        {
            GenePhaseRegion newRegion = newRegions.get(0);
            newRegions.remove(0);

            if(newRegion.length() <= 0)
                continue;

            int index = 0;
            boolean regionSplit = false;
            while(index < regions.size())
            {
                GenePhaseRegion region = regions.get(index);

                if(newRegion == null || region == null)
                {
                    LOGGER.error("isNull: region({}) or new region({})", region == null, newRegion == null);
                    return;
                }

                if (!haveOverlap(region, newRegion, 0))
                {
                    ++index;
                    continue;
                }

                // preserve pre-gene region status only if both regions being combined are pre-gene, otherwise the region will
                // function like a normal phasing region

                // if the region matches on combined phase exactly, then expand the new region to cover both, remove the existing and continue
                if (region.getCombinedPhase() == newRegion.getCombinedPhase() && region.getCombinedPreGeneStatus() == newRegion.getCombinedPreGeneStatus())
                {
                    newRegion.setStart(min(region.start(), newRegion.start()));
                    newRegion.setEnd(max(region.end(), newRegion.end()));
                    newRegion.addPhases(region.getPhaseArray(), region.getPreGenePhaseStatus());

                    regions.remove(index);
                    continue;
                }

                // otherwise split and combine the regions so there are no overlaps
                regionSplit = true;

                // first check for one region enclosing another, the overlaps
                GenePhaseRegion extraRegion = null;

                if (region.start() == newRegion.start() && region.end() == newRegion.end())
                {
                    region.addPhases(newRegion.getPhaseArray(), newRegion.getPreGenePhaseStatus());

                    newRegion.setEnd(newRegion.start()); // won't be added
                }
                else if(region.end() == newRegion.start() || region.start() == newRegion.end())
                {
                    // a precise overlap move base 1 apart and continue
                    if(region.end() == newRegion.start())
                    {
                        region.setEnd(newRegion.start() - 1);
                    }
                    else
                    {
                        newRegion.setEnd(region.start() - 1);
                    }
                }
                else if (region.start() <= newRegion.start() && region.end() >= newRegion.end())
                {
                    // existing region enclosed the new region
                    // split the outer region in 2 and make a combined inner region
                    extraRegion = new GenePhaseRegion(region.GeneId, newRegion.end() + 1, region.end(),
                            region.getPhaseArray(), region.getPreGenePhaseStatus());

                    newRegion.addPhases(region.getPhaseArray(), region.getPreGenePhaseStatus());

                    region.setEnd(newRegion.start() - 1);
                }
                else if (newRegion.start() <= region.start() && newRegion.end() >= region.end())
                {
                    // existing region falls within new region
                    extraRegion = new GenePhaseRegion(newRegion.GeneId, region.end() + 1, newRegion.end(),
                            newRegion.getPhaseArray(), newRegion.getPreGenePhaseStatus());

                    region.addPhases(newRegion.getPhaseArray(), newRegion.getPreGenePhaseStatus());

                    newRegion.setEnd(region.start() - 1);

                }
                else if (newRegion.start() <= region.start())
                {
                    // new region precedes and overlaps the existing

                    // create a new region for the overlapping section
                    extraRegion = new GenePhaseRegion(newRegion.GeneId, region.start(), newRegion.end(),
                            newRegion.getPhaseArray(), newRegion.getPreGenePhaseStatus());

                    extraRegion.addPhases(region.getPhaseArray(), region.getPreGenePhaseStatus());

                    long regionStart = newRegion.end() + 1;
                    long newRegionEnd = region.start() - 1;

                    region.setStart(regionStart);

                    newRegion.setEnd(newRegionEnd);
                }
                else if (region.start() <= newRegion.start())
                {
                    // existing region precedes the new region
                    extraRegion = new GenePhaseRegion(newRegion.GeneId, newRegion.start(), region.end(),
                            newRegion.getPhaseArray(), newRegion.getPreGenePhaseStatus());

                    extraRegion.addPhases(region.getPhaseArray(), region.getPreGenePhaseStatus());

                    long regionEnd = newRegion.start() - 1;
                    long newRegionStart = region.end() + 1;

                    region.setEnd(regionEnd);

                    newRegion.setStart(newRegionStart);
                }

                if(newRegion != null && newRegion.length() >= 1)
                {
                    newRegions.add(newRegion);
                }

                if(extraRegion != null && extraRegion.length() >= 1)
                {
                    newRegions.add(extraRegion);
                }

                ++index;
            }

            if(!regionSplit)
            {
                addRegionInOrder(newRegion, regions);
            }
        }

        // remove 0-1 bases regions
        int index = 0;
        while(index < regions.size())
        {
            if(regions.get(index).length() <= 1)
            {
                regions.remove(index);
            }
            else
            {
                ++index;
            }
        }
    }

    private static void addRegionInOrder(final GenePhaseRegion region, List<GenePhaseRegion> regions)
    {
        int i = 0;
        while(i < regions.size())
        {
            if(region.start() < regions.get(i).start())
                break;

            ++i;
        }

        regions.add(i, region);
    }

    public static void mergePhaseRegions(List<GenePhaseRegion> regions)
    {
        // assumes regions are not overlapping and in order
        int i = 0;
        while(i < regions.size() - 1)
        {
            GenePhaseRegion region1 = regions.get(i);
            GenePhaseRegion region2 = regions.get(i + 1);

            if(region1.end() == region2.start() - 1)
            {
                if (region1.getCombinedPhase() == region2.getCombinedPhase()
                && region1.getCombinedPreGeneStatus() == region2.getCombinedPreGeneStatus())
                {
                    region1.setEnd(region2.end());
                    regions.remove(i + 1);
                    continue;
                }
            }

            ++i;
        }
    }

    public static boolean validateSimpleVsCombinedPhaseRegions(
            final String geneId, final List<GenePhaseRegion> simpleRegions, final List<GenePhaseRegion> combinedRegions)
    {
        boolean hasMismatches = false;

        for(int j = 0; j <= 1; ++j)
        {
            boolean includePreGene = (j == 0);
            long[] phaseCounts = new long[PHASE_MAX];
            long[] combinedPhaseCounts = new long[PHASE_MAX];

            simpleRegions.stream().forEach(x -> x.populateLengthCounts(phaseCounts, includePreGene));
            combinedRegions.stream().forEach(x -> x.populateLengthCounts(combinedPhaseCounts, includePreGene));

            for (int i = 0; i < PHASE_MAX; ++i)
            {
                if (phaseCounts[i] == 0)
                    continue;

                if (abs(phaseCounts[i] - combinedPhaseCounts[i]) / (double) phaseCounts[i] > 0.2)
                {
                    hasMismatches = true;
                    double diffPercent = abs(phaseCounts[i] - combinedPhaseCounts[i]) / (double) phaseCounts[i];

                    LOGGER.warn("geneId({}) phase({}: {}) {} pre-gene simple({}) vs combined({}) diff({})",
                            geneId, i, intAsType(i), includePreGene ? "incl" : "no",
                            phaseCounts[i], combinedPhaseCounts[i], String.format("%.4f", diffPercent));
                }
            }
        }

        if(hasMismatches && combinedRegions.size() <= 10)
        {
            LOGGER.info("geneId({}) simple phases:", geneId);

            for(int i = 0; i < simpleRegions.size(); ++i)
            {
                GenePhaseRegion region = simpleRegions.get(i);
                LOGGER.info("{}: {}", i, region.toString());
            }

            LOGGER.info("geneId({}) combined phases:", geneId);

            for(int i = 0; i < combinedRegions.size(); ++i)
            {
                GenePhaseRegion region = combinedRegions.get(i);
                LOGGER.info("{}: {}", i, region.toString());
            }
        }

        return hasMismatches;
    }
}
