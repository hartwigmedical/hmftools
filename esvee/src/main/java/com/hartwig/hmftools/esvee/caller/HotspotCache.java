package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.sv.LineElements.POLY_A_HOMOLOGY;
import static com.hartwig.hmftools.common.sv.LineElements.POLY_T_HOMOLOGY;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.variant.variantcontext.VariantContext;

public class HotspotCache
{
    private final Map<String,List<KnownHotspot>> mHotspotRegions; // keyed by chromosome start

    private static final String KNOWN_HOTSPOT_FILE = "known_hotspot_file";

    public HotspotCache(final ConfigBuilder configBuilder)
    {
        mHotspotRegions = Maps.newHashMap();

        if(configBuilder.hasValue(KNOWN_HOTSPOT_FILE))
        {
            loadFile(configBuilder.getValue(KNOWN_HOTSPOT_FILE));
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

    public boolean isHotspotVariant(final Variant var)
    {
        // check hotspot rescue
        if(var.isSgl())
            return false;

        if(isPolyATSequence(var.contextStart()) || (var.contextEnd() != null && isPolyATSequence(var.contextEnd())))
            return false;

        return matchesHotspot(var.chromosomeStart(), var.chromosomeEnd(), var.posStart(), var.posEnd(), var.orientStart(), var.orientEnd());
    }

    public boolean isHotspotVariant(final StructuralVariant sv)
    {
        // check hotspot rescue
        if(sv.type() == SGL)
            return false;

        if(isPolyATSequence(sv.startContext()) || (sv.endContext() != null && isPolyATSequence(sv.endContext())))
            return false;

        return matchesHotspot(
                sv.chromosome(true), sv.chromosome(false),
                sv.position(true).intValue(), sv.position(false).intValue(),
                Orientation.fromByte(sv.orientation(true)), Orientation.fromByte(sv.orientation(false)));
    }


    private static boolean isPolyATSequence(final VariantContext variant)
    {
        final String homology = variant.getAttributeAsString(HOMSEQ, "");
        return homology.contains(POLY_A_HOMOLOGY) || homology.contains(POLY_T_HOMOLOGY);
    }

    private boolean matchesHotspot(
            final String chrStart, final String chrEnd, int posStart, int posEnd, Orientation orientStart, Orientation orientEnd)
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

                // note BED file adjustment for start position
                ChrBaseRegion regionStart = new ChrBaseRegion(chrStart, Integer.parseInt(items[1]) + 1, Integer.parseInt(items[2]));
                ChrBaseRegion regionEnd = new ChrBaseRegion(chrEnd, Integer.parseInt(items[4]) + 1, Integer.parseInt(items[5]));
                Orientation orientStart = Orientation.fromChar(items[8].charAt(0));
                Orientation orientEnd = Orientation.fromChar(items[9].charAt(0));
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

            SV_LOGGER.info("loaded {} known hotspot records from file", itemCount, filename);
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to load hotspot data file({}): {}", filename, e.toString());
        }
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(KNOWN_HOTSPOT_FILE, false, "Known fusion BED file");
    }
}
