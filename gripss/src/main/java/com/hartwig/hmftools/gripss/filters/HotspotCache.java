package com.hartwig.hmftools.gripss.filters;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssConfig.GR_LOGGER;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_HOMSEQ;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.POLY_A;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.POLY_T;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.links.LinkRescue;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

import htsjdk.variant.variantcontext.VariantContext;

public class HotspotCache
{
    private final Map<String,List<KnownHotspot>> mHotspotRegions; // keyed by chromosome start

    private static final String KNOWN_HOTSPOT_FILE = "known_hotspot_file";

    public HotspotCache(final CommandLine cmd)
    {
        mHotspotRegions = Maps.newHashMap();

        if(cmd != null)
        {
            loadFile(cmd.getOptionValue(KNOWN_HOTSPOT_FILE));
        }
    }

    public void addHotspot(final KnownHotspot hotspot)
    {
        List<KnownHotspot> hotspots = mHotspotRegions.get(hotspot.RegionStart.Chromosome);

        if(hotspots == null)
        {
            hotspots = Lists.newArrayList();
            mHotspotRegions.put(hotspot.RegionStart.Chromosome, hotspots);
        }

        hotspots.add(hotspot);
    }

    public boolean isHotspotVariant(final SvData sv)
    {
        // check hotspot rescue
        if(sv.isSgl())
            return false;

        if(LinkRescue.tooShortToRescue(sv.type(), sv.length()))
            return false;

        if(isPolyATSequence(sv.contextStart()) || (sv.contextEnd() != null && isPolyATSequence(sv.contextEnd())))
            return false;

        return matchesHotspot(sv.chromosomeStart(), sv.chromosomeEnd(), sv.posStart(), sv.posEnd(), sv.orientStart(), sv.orientEnd());
    }

    public boolean isHotspotVariant(final StructuralVariant sv)
    {
        // check hotspot rescue
        if(sv.type() == SGL)
            return false;

        if(LinkRescue.tooShortToRescue(sv.type(), SvData.length(sv)))
            return false;

        else if(isPolyATSequence(sv.startContext()) || (sv.endContext() != null && isPolyATSequence(sv.endContext())))
            return false;

        return matchesHotspot(
                sv.chromosome(true), sv.chromosome(false),
                sv.position(true).intValue(), sv.position(false).intValue(),
                sv.orientation(true), sv.orientation(false));
    }


    private static boolean isPolyATSequence(final VariantContext variant)
    {
        final String homology = variant.getAttributeAsString(VT_HOMSEQ, "");
        return homology.contains(POLY_A) || homology.contains(POLY_T);
    }

    private boolean matchesHotspot(final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd)
    {
        List<KnownHotspot> regions = mHotspotRegions.get(chrStart);
        if(regions != null)
        {
            for(KnownHotspot region : regions)
            {
                if(region.matches(chrStart, chrEnd, posStart, posEnd, orientStart, orientEnd))
                    return true;
            }
        }

        return false;
    }

    public boolean matchesHotspotBreakend(final String chromosome, int position)
    {
        List<KnownHotspot> regions = mHotspotRegions.get(chromosome);

        if(regions != null)
        {
            for(KnownHotspot region : regions)
            {
                if(region.RegionStart.containsPosition(position) || region.RegionEnd.containsPosition(position))
                    return true;
            }
        }

        return false;
    }

    private void loadFile(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            int itemCount = 0;
            String line = fileReader.readLine();

            while(line != null)
            {
                final String[] items = line.split("\t", -1);

                String chrStart = items[0];
                String chrEnd = items[3];

                ChrBaseRegion regionStart = new ChrBaseRegion(chrStart, Integer.parseInt(items[1]), Integer.parseInt(items[2]));
                ChrBaseRegion regionEnd = new ChrBaseRegion(chrEnd, Integer.parseInt(items[4]), Integer.parseInt(items[5]));
                Byte orientStart = items[8].equals("+") ? POS_ORIENT : NEG_ORIENT;
                Byte orientEnd = items[9].equals("+") ? POS_ORIENT : NEG_ORIENT;
                String geneInfo = items[6];

                KnownHotspot knownHotspot = new KnownHotspot(regionStart, orientStart, regionEnd, orientEnd, geneInfo);

                // add to both chromosome lists
                List<KnownHotspot> svRegions = mHotspotRegions.get(chrStart);

                if(svRegions == null)
                {
                    svRegions = Lists.newArrayList();
                    mHotspotRegions.put(chrStart, svRegions);
                }

                svRegions.add(knownHotspot);

                if(!chrStart.equals(chrEnd))
                {
                    svRegions = mHotspotRegions.get(chrEnd);

                    if(svRegions == null)
                    {
                        svRegions = Lists.newArrayList();
                        mHotspotRegions.put(chrEnd, svRegions);
                    }

                    svRegions.add(knownHotspot);
                }

                ++itemCount;

                line = fileReader.readLine();
            }

            GR_LOGGER.info("loaded {} known hotspot records from file", itemCount, filename);
        }
        catch(IOException e)
        {
            GR_LOGGER.error("failed to load hotspot data file({}): {}", filename, e.toString());
            return;
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(KNOWN_HOTSPOT_FILE, true, "Known fusion BED file");
    }
}
