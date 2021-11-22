package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssConfig.GR_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.gripss.common.SvData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class PonCache
{
    private final Map<String,List<PonSvRegion>> mSvRegions;
    private final Map<String,List<PonSglRegion>> mSglRegions;
    private final List<String> mRestrictedChromosomes;
    private final int mPositionMargin;

    private static final String GERMLINE_PON_BED_SV_FILE = "pon_sv_file";
    private static final String GERMLINE_PON_BED_SGL_FILE = "pon_sgl_file";
    private static final String GERMLINE_PON_MARGIN = "pon_margin";

    public PonCache(final CommandLine cmd, final List<String> restrictedChromosomes)
    {
        mSvRegions = Maps.newHashMap();
        mSglRegions = Maps.newHashMap();
        mRestrictedChromosomes = restrictedChromosomes;

        if(cmd != null)
        {
            mPositionMargin = Integer.parseInt(cmd.getOptionValue(GERMLINE_PON_MARGIN, "0"));
            loadPonSvFile(cmd.getOptionValue(GERMLINE_PON_BED_SV_FILE));
            loadPonSglFile(cmd.getOptionValue(GERMLINE_PON_BED_SGL_FILE));
        }
        else
        {
            mPositionMargin = 0;
        }
    }

    public int getPonCount(final SvData var)
    {
        if(var.isSgl())
        {
            List<PonSglRegion> regions = mSglRegions.get(var.chromosomeStart());
            if(regions != null)
            {
                return findSglPonMatch(regions, var);
            }
        }
        else
        {
            List<PonSvRegion> regions = mSvRegions.get(var.chromosomeStart());
            if(regions != null)
            {
                return findPonMatch(regions, var);
            }
        }

        return 0;
    }

    private int findPonMatch(final List<PonSvRegion> regions, final SvData var)
    {
        // use a binary search
        int posStart = var.posStart();
        int currentIndex = regions.size() / 2;
        int lowerIndex = 0;
        int upperIndex = regions.size() - 1;
        while(true)
        {
            PonSvRegion region = regions.get(currentIndex);

            if(region.withinStartRegion(posStart, mPositionMargin))
            {
                // test the PON entries around this position
                return findPonMatch(regions, var, currentIndex);
            }

            if(region.RegionStart.end() < posStart)
            {
                // SV's position is higher than the current index
                lowerIndex = currentIndex;
            }
            else
            {
                upperIndex = currentIndex;
            }

            int nextIndex = (lowerIndex + upperIndex) / 2;
            if(currentIndex == nextIndex)
                return 0;

            currentIndex = nextIndex;
        }
    }

    private int findPonMatch(final List<PonSvRegion> regions, final SvData var, int startIndex)
    {
        // search and up and down from this entry point for a PON match
        int posStart = var.posStart();

        for(int i = 0; i <= 1; ++i)
        {
            boolean searchUp = (i == 0);
            int currentIndex = startIndex;

            while(currentIndex >= 0 && currentIndex < regions.size())
            {
                PonSvRegion region = regions.get(currentIndex);

                if(!region.withinStartRegion(posStart, mPositionMargin))
                    break;

                if(region.matches(var, mPositionMargin))
                    return region.PonCount;

                if(searchUp)
                    ++currentIndex;
                else
                    --currentIndex;
            }
        }

        return 0;
    }

    private int findSglPonMatch(final List<PonSglRegion> regions, final SvData var)
    {
        // use a binary search
        int posStart = var.posStart();
        int currentIndex = regions.size() / 2;
        int lowerIndex = 0;
        int upperIndex = regions.size() - 1;
        while(true)
        {
            PonSglRegion region = regions.get(currentIndex);

            if(region.withinRegion(posStart, mPositionMargin))
            {
                // test the PON entries around this position
                return findSglPonMatch(regions, var, currentIndex);
            }

            if(region.Region.start() < posStart)
            {
                // SV's position is higher than the current index
                lowerIndex = currentIndex;
            }
            else
            {
                upperIndex = currentIndex;
            }

            int nextIndex = (lowerIndex + upperIndex) / 2;
            if(currentIndex == nextIndex)
                return 0;

            currentIndex = nextIndex;
        }
    }

    private int findSglPonMatch(final List<PonSglRegion> regions, final SvData var, int startIndex)
    {
        // search and up and down from this entry point for a PON match
        int posStart = var.posStart();

        for(int i = 0; i <= 1; ++i)
        {
            boolean searchUp = (i == 0);
            int currentIndex = startIndex;

            while(currentIndex >= 0 && currentIndex < regions.size())
            {
                PonSglRegion region = regions.get(currentIndex);

                if(!region.withinRegion(posStart, mPositionMargin))
                    break;

                if(region.matches(var, mPositionMargin))
                    return region.PonCount;

                if(searchUp)
                    ++currentIndex;
                else
                    --currentIndex;
            }
        }

        return 0;
    }

    private void loadPonSvFile(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            int itemCount = 0;
            String line = null;
            String currentChr = "";
            List<PonSvRegion> svRegions = null;
            ChrBaseRegion lastRegion = null;

            while((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split("\t", -1);

                String chrStart = items[0];
                String chrEnd = items[3];

                if(!mRestrictedChromosomes.isEmpty())
                {
                    if(!mRestrictedChromosomes.contains(chrStart) || !mRestrictedChromosomes.contains(chrEnd))
                        continue;
                }

                if(!chrStart.equals(currentChr))
                {
                    currentChr = chrStart;
                    svRegions = Lists.newArrayList();
                    mSvRegions.put(chrStart, svRegions);
                    lastRegion = null;
                }

                ChrBaseRegion regionStart = new ChrBaseRegion(chrStart, Integer.parseInt(items[1]), Integer.parseInt(items[2]));
                ChrBaseRegion regionEnd = new ChrBaseRegion(chrEnd, Integer.parseInt(items[4]), Integer.parseInt(items[5]));

                Byte orientStart = items[8].equals("+") ? POS_ORIENT : NEG_ORIENT;
                Byte orientEnd = items[9].equals("+") ? POS_ORIENT : NEG_ORIENT;
                int ponCount = Integer.parseInt(items[7]);

                svRegions.add(new PonSvRegion(regionStart, orientStart, regionEnd, orientEnd, ponCount));
                ++itemCount;

                if(lastRegion != null && lastRegion.end() > regionStart.end())
                {
                    GR_LOGGER.warn("SV PON not ordered: last({}) vs this({})", lastRegion, regionStart);
                }

                lastRegion = regionStart;
            }

            GR_LOGGER.info("loaded {} germline SV PON records from file({})", itemCount, filename);
        }
        catch(IOException e)
        {
            GR_LOGGER.error("failed to load germline SV PON file({}): {}", filename, e.toString());
            return;
        }
    }

    private void loadPonSglFile(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            int itemCount = 0;
            String line = null;
            String currentChr = "";
            List<PonSglRegion> sglRegions = null;
            BaseRegion lastRegion = null;

            while((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split("\t", -1);

                String chr = items[0];

                if(!mRestrictedChromosomes.isEmpty() && !mRestrictedChromosomes.contains(chr))
                    continue;

                if(!chr.equals(currentChr))
                {
                    currentChr = chr;
                    sglRegions = Lists.newArrayList();
                    mSglRegions.put(chr, sglRegions);
                    lastRegion = null;
                }

                BaseRegion region = new BaseRegion(Integer.parseInt(items[1]), Integer.parseInt(items[2]));

                Byte orient = items[5].equals("+") ? POS_ORIENT : NEG_ORIENT;
                int ponCount = Integer.parseInt(items[4]);

                sglRegions.add(new PonSglRegion(region, orient, ponCount));
                ++itemCount;
                
                if(lastRegion != null && lastRegion.start() > region.start())
                {
                    GR_LOGGER.warn("SGL PON not ordered: last({}) vs this({})", lastRegion, region);
                }

                lastRegion = region;
            }

            GR_LOGGER.info("loaded {} germline SGL PON records from file", itemCount, filename);
        }
        catch(IOException e)
        {
            GR_LOGGER.error("failed to load germline SGL PON file({}): {}", filename, e.toString());
            return;
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(GERMLINE_PON_BED_SV_FILE, true, "PON for SV positions");
        options.addOption(GERMLINE_PON_BED_SGL_FILE, true, "PON for SGL positions");
        options.addOption(GERMLINE_PON_MARGIN, true, "PON permitted matching position margin");
    }

    private class PonSvRegion
    {
        public final ChrBaseRegion RegionStart;
        public final Byte OrientStart;
        public final ChrBaseRegion RegionEnd;
        public final Byte OrientEnd;
        public final int PonCount;

        public PonSvRegion(
                final ChrBaseRegion regionStart, final Byte orientStart, final ChrBaseRegion regionEnd, final Byte orientEnd, final int ponCount)
        {
            RegionStart = regionStart;
            OrientStart = orientStart;
            RegionEnd = regionEnd;
            OrientEnd = orientEnd;
            PonCount = ponCount;
        }

        public boolean withinStartRegion(int position, int margin)
        {
            return position >= RegionStart.start() - margin && position <= RegionStart.end() + margin;
        }

        public boolean matches(final SvData var, int margin)
        {
            if(var.posStart() < RegionStart.start() - margin || var.posStart() > RegionStart.end() + margin)
                return false;

            if(var.posEnd() < RegionEnd.start() - margin || var.posEnd() > RegionEnd.end() + margin)
                return false;

            return OrientStart == var.orientStart() && OrientEnd == var.orientEnd();
        }
    }

    private class PonSglRegion
    {
        public final BaseRegion Region;
        public final Byte Orient;
        public final int PonCount;

        public PonSglRegion(final BaseRegion region, final Byte orient, final int ponCount)
        {
            Region = region;
            Orient = orient;
            PonCount = ponCount;
        }

        public boolean withinRegion(int position, int margin)
        {
            return position >= Region.start() - margin && position <= Region.end() + margin;
        }

        public boolean matches(final SvData var, int margin)
        {
            return withinRegion(var.posStart(), margin) && Orient == var.orientStart();
        }

    }
}
