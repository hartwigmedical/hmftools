package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class HotspotCache
{
    private final Map<String, List<KnownHotspot>> mHotspotRegions; // keyed by chromosome start

    public HotspotCache(final String filename)
    {
        mHotspotRegions = Maps.newHashMap();
        loadFile(filename);
    }

    public boolean matchesHotspot(final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd)
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

            SV_LOGGER.info("loaded {} known hotspot records from file", itemCount, filename);
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to load hotspot data file({}): {}", filename, e.toString());
            return;
        }
    }

    public class KnownHotspot
    {
        public final ChrBaseRegion RegionStart;
        public final Byte OrientStart;
        public final ChrBaseRegion RegionEnd;
        public final Byte OrientEnd;
        public final String GeneInfo;

        public KnownHotspot(
                final ChrBaseRegion regionStart, final Byte orientStart, final ChrBaseRegion regionEnd, final Byte orientEnd,
                final String geneInfo)
        {
            RegionStart = regionStart;
            OrientStart = orientStart;
            RegionEnd = regionEnd;
            OrientEnd = orientEnd;
            GeneInfo = geneInfo;
        }

        public boolean matches(final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd)
        {
            if(RegionStart.containsPosition(chrStart, posStart) && OrientStart == orientStart
            && RegionEnd.containsPosition(chrEnd, posEnd) && OrientEnd == orientEnd)
            {
                return true;
            }

            if(RegionStart.containsPosition(chrEnd, posEnd) && OrientStart == orientEnd
            && RegionEnd.containsPosition(chrStart, posStart) && OrientEnd == orientStart)
            {
                return true;
            }

            return false;
        }
    }
}
