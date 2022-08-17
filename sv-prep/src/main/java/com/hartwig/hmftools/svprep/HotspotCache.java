package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class HotspotCache
{
    private final Map<String,List<KnownHotspot>> mHotspotRegions; // keyed by chromosome start
    private final boolean mIsValid;

    public HotspotCache(final String filename)
    {
        mHotspotRegions = Maps.newHashMap();
        mIsValid = loadFile(filename);
    }

    public boolean isValid() { return mIsValid; }

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

    public List<ChrBaseRegion> findMatchingRegions(final ChrBaseRegion region)
    {
        List<ChrBaseRegion> matchedRegions = Lists.newArrayList();

        List<KnownHotspot> knownHotspots = mHotspotRegions.get(region.Chromosome);

        if(knownHotspots != null)
        {
            for(KnownHotspot knownHotspot : knownHotspots)
            {
                if(knownHotspot.RegionStart.overlaps(region))
                    matchedRegions.add(knownHotspot.RegionStart);

                if(knownHotspot.RegionEnd.overlaps(region))
                    matchedRegions.add(knownHotspot.RegionEnd);
            }
        }

        return matchedRegions;
    }

    private boolean loadFile(final String filename)
    {
        if(filename == null)
            return true;

        try
        {
            BufferedReader fileReader = createBufferedReader(filename);

            int itemCount = 0;
            String line = "";

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split("\t", -1);

                if(values.length < 10)
                {
                    SV_LOGGER.error("invalid hotspot entry: {}", line);
                    return false;
                }

                String chrStart = values[0];
                String chrEnd = values[3];

                ChrBaseRegion regionStart = new ChrBaseRegion(chrStart, Integer.parseInt(values[1]) + 1, Integer.parseInt(values[2]));
                ChrBaseRegion regionEnd = new ChrBaseRegion(chrEnd, Integer.parseInt(values[4]) + 1, Integer.parseInt(values[5]));
                Byte orientStart = values[8].equals("+") ? POS_ORIENT : NEG_ORIENT;
                Byte orientEnd = values[9].equals("+") ? POS_ORIENT : NEG_ORIENT;
                String geneInfo = values[6];

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
            }

            SV_LOGGER.info("loaded {} known hotspot records from file", itemCount, filename);
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to load hotspot data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
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
