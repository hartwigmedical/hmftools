package com.hartwig.hmftools.linx.germline;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class GermlinePonCache
{
    private final Map<String,List<PonSvRegion>> mSvRegions;
    private final Map<String,List<PonSglRegion>> mSglRegions;

    private static final String GERMLINE_PON_BED_SV_FILE = "germine_pon_sv_file";
    private static final String GERMLINE_PON_BED_SGL_FILE = "germine_pon_sgl_file";

    public GermlinePonCache(final CommandLine cmd)
    {
        mSvRegions = Maps.newHashMap();
        mSglRegions = Maps.newHashMap();

        if(cmd != null)
        {
            loadPonSvFile(cmd.getOptionValue(GERMLINE_PON_BED_SV_FILE));
            loadPonSglFile(cmd.getOptionValue(GERMLINE_PON_BED_SGL_FILE));
        }
    }

    public int getPonCount(final SvVarData var)
    {
        if(var.isSglBreakend())
        {
            List<PonSglRegion> regions = mSglRegions.get(var.chromosome(true));
            if(regions != null)
            {
                for(PonSglRegion region : regions)
                {
                    if(region.matches(var))
                        return region.PonCount;
                }
            }
        }
        else
        {
            List<PonSvRegion> regions = mSvRegions.get(var.chromosome(true));
            if(regions != null)
            {
                for(PonSvRegion region : regions)
                {
                    if(region.matches(var))
                        return region.PonCount;
                }
            }
        }

        return 0;
    }

    private void loadPonSvFile(final String filename)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            int itemCount = 0;
            String line = fileReader.readLine();
            String currentChr = "";
            List<PonSvRegion> svRegions = null;

            while(line != null)
            {
                final String[] items = line.split("\t", -1);

                String chrStart = items[0];

                if(!chrStart.equals(currentChr))
                {
                    currentChr = chrStart;
                    svRegions = Lists.newArrayList();
                    mSvRegions.put(chrStart, svRegions);
                }

                ChrBaseRegion regionStart = new ChrBaseRegion(chrStart, Integer.parseInt(items[1]), Integer.parseInt(items[2]));
                ChrBaseRegion regionEnd = new ChrBaseRegion(items[3], Integer.parseInt(items[4]), Integer.parseInt(items[5]));
                Byte orientStart = items[8].equals("+") ? POS_ORIENT : NEG_ORIENT;
                Byte orientEnd = items[9].equals("+") ? POS_ORIENT : NEG_ORIENT;
                int ponCount = Integer.parseInt(items[7]);

                svRegions.add(new PonSvRegion(regionStart, orientStart, regionEnd, orientEnd, ponCount));
                ++itemCount;

                line = fileReader.readLine();
            }

            LNX_LOGGER.info("loaded {} germline SV PON records from file", itemCount, filename);
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to load germline SV PON file({}): {}", filename, e.toString());
            return;
        }
    }

    private void loadPonSglFile(final String filename)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            int itemCount = 0;
            String line = fileReader.readLine();
            String currentChr = "";
            List<PonSglRegion> sglRegions = null;

            while(line != null)
            {
                final String[] items = line.split("\t", -1);

                String chr = items[0];

                if(!chr.equals(currentChr))
                {
                    currentChr = chr;
                    sglRegions = Lists.newArrayList();
                    mSglRegions.put(chr, sglRegions);
                }

                BaseRegion region = new BaseRegion(Integer.parseInt(items[1]), Integer.parseInt(items[2]));
                Byte orient = items[5].equals("+") ? POS_ORIENT : NEG_ORIENT;
                int ponCount = Integer.parseInt(items[4]);

                sglRegions.add(new PonSglRegion(region, orient, ponCount));
                ++itemCount;

                line = fileReader.readLine();
            }

            LNX_LOGGER.info("loaded {} germline SGL PON records from file", itemCount, filename);
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to load germline SGL PON file({}): {}", filename, e.toString());
            return;
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(GERMLINE_PON_BED_SV_FILE, true, "Optional: write all batch-run output files");
        options.addOption(GERMLINE_PON_BED_SGL_FILE, true, "Optional: include all SV table fields (batch-mode)");
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

        public boolean matches(final SvVarData var)
        {
            return RegionStart.containsPosition(var.chromosome(true), var.position(true))
                    && OrientStart == var.orientation(true)
                    && RegionEnd.containsPosition(var.chromosome(false), var.position(false))
                    && OrientEnd == var.orientation(false);
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

        public boolean matches(final SvVarData var)
        {
            return Region.containsPosition(var.position(true)) && Orient == var.orientation(true);
        }

    }
}
