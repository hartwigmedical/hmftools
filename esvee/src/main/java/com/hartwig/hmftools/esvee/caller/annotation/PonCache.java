package com.hartwig.hmftools.esvee.caller.annotation;

import static java.lang.Integer.max;
import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.DEFAULT_PON_DISTANCE;
import static com.hartwig.hmftools.esvee.common.FilterType.PON;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.caller.Breakend;
import com.hartwig.hmftools.esvee.caller.Variant;

public class PonCache
{
    private final Map<String,List<PonSvRegion>> mSvRegions;
    private final Map<String,List<PonSglRegion>> mSglRegions;
    private final int mPositionMargin;
    private final boolean mAllowUnordered;

    // keep indices into the 2 collections assuming that requests to match on the PON will be made sequentially through the genome
    private String mCurrentSvChromosome;
    private int mCurrentSvIndex;
    private String mCurrentSglChromosome;
    private int mCurrentSglIndex;
    private boolean mHasValidData;

    private static final String GERMLINE_PON_BED_SV_FILE = "pon_sv_file";
    private static final String GERMLINE_PON_BED_SGL_FILE = "pon_sgl_file";
    private static final String GERMLINE_PON_MARGIN = "pon_margin";

    public PonCache(final ConfigBuilder configBuilder)
    {
        this(configBuilder.getInteger(GERMLINE_PON_MARGIN), configBuilder.getValue(GERMLINE_PON_BED_SV_FILE),
                configBuilder.getValue(GERMLINE_PON_BED_SGL_FILE), false);
    }

    public PonCache(final int margin, final String ponSvFile, final String ponSglFile, boolean allowUnordered)
    {
        mSvRegions = Maps.newHashMap();
        mSglRegions = Maps.newHashMap();
        mAllowUnordered = allowUnordered;
        mHasValidData = true;

        mPositionMargin = margin;

        if(ponSvFile != null)
            loadPonSvFile(ponSvFile);

        if(ponSglFile != null)
            loadPonSglFile(ponSglFile);

        mCurrentSglChromosome = "";
        mCurrentSvChromosome = "";
        mCurrentSglIndex = 0;
        mCurrentSvIndex = 0;
    }

    public boolean hasValidData() { return mHasValidData; }

    public Map<String,List<PonSvRegion>> svRegions() { return mSvRegions; }
    public Map<String,List<PonSglRegion>> sglRegions() { return mSglRegions; }

    public void annotateVariants(final List<Variant> variantList)
    {
        for(Variant var : variantList)
        {
            int ponCount = getPonCount(var);

            if(ponCount > 0)
            {
                var.setPonCount(ponCount);

                var.addFilter(PON);
            }
        }
    }

    public int getPonCount(final Variant var)
    {
        // matching routine:
        // - get regions by chromosome
        // - use a binary search using start position
        // - if a region is matched, search up and down from their checking both positions and orientations

        if(var.isSgl())
        {
            List<PonSglRegion> regions = mSglRegions.get(var.chromosomeStart());
            if(regions != null)
            {
                if(!mCurrentSglChromosome.equals(var.chromosomeStart()))
                {
                    mCurrentSglChromosome = var.chromosomeStart();
                    mCurrentSglIndex = 0;
                }

                return findSglPonMatch(regions, var);
            }
        }
        else
        {
            List<PonSvRegion> regions = mSvRegions.get(var.chromosomeStart());
            if(regions != null)
            {
                if(!mCurrentSvChromosome.equals(var.chromosomeStart()))
                {
                    mCurrentSvChromosome = var.chromosomeStart();
                    mCurrentSvIndex = 0;
                }

                return findPonMatch(regions, var);
            }
        }

        return 0;
    }

    private static int[] breakendMargin(final Breakend breakend)
    {
        int[] margins = new int[SE_PAIR];
        int inexactHomology = abs(breakend.IsStart ? breakend.InexactHomology.Start : breakend.InexactHomology.End);
        margins[SE_START] = min(breakend.ConfidenceInterval.Start, -inexactHomology);
        margins[SE_END] = max(breakend.ConfidenceInterval.End, inexactHomology);
        return margins;
    }

    private int findPonMatch(final List<PonSvRegion> regions, final Variant var)
    {
        final int[] marginStart = breakendMargin(var.breakendStart());
        final int[] marginEnd = breakendMargin(var.breakendEnd());

        BaseRegion svStart = new BaseRegion(
                var.posStart() + marginStart[SE_START] - mPositionMargin,
                var.posStart() + marginStart[SE_END] + mPositionMargin);

        for(; mCurrentSvIndex < regions.size(); ++mCurrentSvIndex)
        {
            PonSvRegion region = regions.get(mCurrentSvIndex);

            if(region.RegionStart.overlaps(svStart))
            {
                // test the PON entries around this position
                ChrBaseRegion svEnd = new ChrBaseRegion(
                        var.chromosomeEnd(),
                        var.posEnd() + marginEnd[SE_START] - mPositionMargin,
                        var.posEnd() + marginEnd[SE_END] + mPositionMargin);

                return findPonMatch(regions, var, svStart, svEnd, mCurrentSvIndex);
            }

            // exit if the PON is now past this point and retreat one position
            if(region.RegionStart.start() > svStart.end())
                break;
        }

        if(mCurrentSvIndex > 0)
            --mCurrentSvIndex;

        return 0;
    }

    private int findPonMatch(
            final List<PonSvRegion> regions, final Variant var, final BaseRegion svStart, final ChrBaseRegion svEnd, int startIndex)
    {
        // search and up and down from this entry point for a PON match
        for(int i = 0; i <= 1; ++i)
        {
            boolean searchUp = (i == 0);
            int currentIndex = startIndex;

            while(currentIndex >= 0 && currentIndex < regions.size())
            {
                PonSvRegion region = regions.get(currentIndex);

                if(searchUp && region.RegionStart.start() > svStart.end())
                    break;

                if(!searchUp && region.RegionStart.end() < svStart.start())
                    break;

                if(region.matches(svStart, svEnd, var.orientStart(), var.orientEnd()))
                    return region.PonCount;

                if(searchUp)
                    ++currentIndex;
                else
                    --currentIndex;
            }
        }

        return 0;
    }

    private int findSglPonMatch(final List<PonSglRegion> regions, final Variant var)
    {
        final int[] marginStart = breakendMargin(var.breakendStart());

        BaseRegion svStart = new BaseRegion(
                var.posStart() + marginStart[SE_START] - mPositionMargin,
                var.posStart() + marginStart[SE_END] + mPositionMargin);

        for(; mCurrentSglIndex < regions.size(); ++mCurrentSglIndex)
        {
            PonSglRegion region = regions.get(mCurrentSglIndex);

            if(region.Region.overlaps(svStart))
            {
                // test the PON entries around this position
                return findSglPonMatch(regions, var, svStart, mCurrentSglIndex);
            }

            // exit if the PON is now past this point and retreat one position
            if(region.Region.start() > svStart.end())
            {
                if(mCurrentSglIndex > 0)
                    --mCurrentSglIndex;

                break;
            }
        }

        return 0;
    }

    private int findSglPonMatch(final List<PonSglRegion> regions, final Variant var, final BaseRegion svStart, int startIndex)
    {
        // search and up and down from this entry point for a PON match

        for(int i = 0; i <= 1; ++i)
        {
            boolean searchUp = (i == 0);
            int currentIndex = startIndex;

            while(currentIndex >= 0 && currentIndex < regions.size())
            {
                PonSglRegion region = regions.get(currentIndex);

                if(searchUp && region.Region.start() > svStart.end())
                    break;

                if(!searchUp && region.Region.end() < svStart.start())
                    break;

                if(region.matches(svStart, var.orientStart()))
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
            BufferedReader fileReader = createBufferedReader(filename);

            int itemCount = 0;
            String line = null;
            String currentChr = "";
            List<PonSvRegion> svRegions = null;
            BaseRegion lastRegion = null;

            // fields: ChrStart,PosStartBegin,PosStartEnd,ChrEnd,PosEndBegin,PosEndEnd,Unknown,PonCount,OrientStart,OrientEnd

            while((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(TSV_DELIM, -1);

                String chrStart = items[0];
                String chrEnd = items[3];

                if(!chrStart.equals(currentChr))
                {
                    currentChr = chrStart;
                    svRegions = Lists.newArrayList();
                    mSvRegions.put(chrStart, svRegions);
                    lastRegion = null;
                }

                // note BED start position adjustment
                BaseRegion regionStart = new BaseRegion(Integer.parseInt(items[1]) + 1, Integer.parseInt(items[2]));
                ChrBaseRegion regionEnd = new ChrBaseRegion(chrEnd, Integer.parseInt(items[4]) + 1, Integer.parseInt(items[5]));

                Orientation orientStart = Orientation.fromChar(items[8].charAt(0));
                Orientation orientEnd = Orientation.fromChar(items[9].charAt(0));
                int ponCount = Integer.parseInt(items[7]);

                svRegions.add(new PonSvRegion(regionStart, orientStart, regionEnd, orientEnd, ponCount));
                ++itemCount;

                if(!mAllowUnordered && lastRegion != null && lastRegion.start() > regionStart.start())
                {
                    SV_LOGGER.warn("SV PON not ordered: last({}) vs this({})", lastRegion, regionStart);
                }

                lastRegion = regionStart;
            }

            SV_LOGGER.info("loaded {} germline SV PON records from file({})", itemCount, filename);
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to load germline SV PON file({}): {}", filename, e.toString());
            mHasValidData = false;
            return;
        }
    }

    private void loadPonSglFile(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            BufferedReader fileReader = createBufferedReader(filename);

            int itemCount = 0;
            String line = null;
            String currentChr = "";
            List<PonSglRegion> sglRegions = null;
            BaseRegion lastRegion = null;

            // fields: Chr,PosBegin,PosEnd,Unknown,PonCount,Orientation

            while((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(TSV_DELIM, -1);

                String chr = items[0];

                if(!chr.equals(currentChr))
                {
                    currentChr = chr;
                    sglRegions = Lists.newArrayList();
                    mSglRegions.put(chr, sglRegions);
                    lastRegion = null;
                }

                BaseRegion region = new BaseRegion(Integer.parseInt(items[1]) + 1, Integer.parseInt(items[2]));

                Orientation orient = Orientation.fromChar(items[5].charAt(0));
                int ponCount = Integer.parseInt(items[4]);

                sglRegions.add(new PonSglRegion(region, orient, ponCount));
                ++itemCount;
                
                if(!mAllowUnordered && lastRegion != null && lastRegion.start() > region.start())
                {
                    SV_LOGGER.warn("SGL PON not ordered: last({}) vs this({})", lastRegion, region);
                }

                lastRegion = region;
            }

            SV_LOGGER.info("loaded {} germline SGL PON records from file({})", itemCount, filename);
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to load germline SGL PON file({}): {}", filename, e.toString());
            mHasValidData = false;
            return;
        }
    }

    public void addPonSvRegion(
            final String chrStart, final BaseRegion regionStart, final Orientation orientStart,
            final ChrBaseRegion regionEnd, final Orientation orientEnd, final int ponCount)
    {
        List<PonSvRegion> regions = mSvRegions.get(chrStart);

        if(regions == null)
        {
            regions = Lists.newArrayList();
            mSvRegions.put(chrStart, regions);
        }

        regions.add(new PonSvRegion(regionStart, orientStart, regionEnd, orientEnd, ponCount));
    }

    public void addPonSglRegion(final String chromosome, BaseRegion region, final Orientation orient, final int ponCount)
    {
        List<PonSglRegion> regions = mSglRegions.get(chromosome);

        if(regions == null)
        {
            regions = Lists.newArrayList();
            mSglRegions.put(chromosome, regions);
        }

        regions.add(new PonSglRegion(region, orient, ponCount));
    }

    public void clear()
    {
        mCurrentSvIndex = 0;
        mCurrentSglIndex = 0;
        mSvRegions.clear();
        mSglRegions.clear();
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(GERMLINE_PON_BED_SV_FILE, false, "PON for SV positions");
        configBuilder.addPath(GERMLINE_PON_BED_SGL_FILE, false, "PON for SGL positions");
        configBuilder.addInteger(
                GERMLINE_PON_MARGIN, "PON permitted matching position margin", DEFAULT_PON_DISTANCE);
    }

}
