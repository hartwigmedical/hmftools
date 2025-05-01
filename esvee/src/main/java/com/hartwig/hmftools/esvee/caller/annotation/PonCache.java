package com.hartwig.hmftools.esvee.caller.annotation;

import static java.lang.Integer.max;
import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sv.SvVcfTags.INSALN;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.MULTI_MAPPED_ALT_ALIGNMENT_REGIONS_V37;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.MULTI_MAPPED_ALT_ALIGNMENT_REGIONS_V38;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.SSX2_REGIONS_V37;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.SSX2_REGIONS_V38;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.DEFAULT_PON_DISTANCE;
import static com.hartwig.hmftools.esvee.common.FilterType.PON;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.assembly.alignment.AlternativeAlignment;
import com.hartwig.hmftools.esvee.caller.Breakend;
import com.hartwig.hmftools.esvee.caller.Variant;

public class PonCache
{
    private final Map<String,List<PonSvRegion>> mSvRegions; // SV entries keyed by starting/lower chromosome
    private final Map<String,List<PonSglRegion>> mSglRegions; // SGL entries keyed by chromosome
    private final int mPositionMargin;
    private final boolean mAllowUnordered;

    // keep indices into the 2 collections assuming that requests to match on the PON will be made sequentially through the genome
    private int mCurrentSvIndex;
    private int mCurrentSglIndex;
    private boolean mHasValidData;

    private final boolean mRunChecks = false;

    private static final String GERMLINE_PON_BED_SV_FILE = "pon_sv_file";
    private static final String GERMLINE_PON_BED_SGL_FILE = "pon_sgl_file";
    public static final String GERMLINE_PON_MARGIN = "pon_margin";

    public static final String ARTEFACT_PON_BED_SV_FILE = "artefact_pon_sv_file";
    public static final String ARTEFACT_PON_BED_SGL_FILE = "artefact_pon_sgl_file";

    public static final String FLD_PON_COUNT = "PonCount";

    private final List<ChrBaseRegion> mSpecificSglFusionRegions;

    public PonCache(final ConfigBuilder configBuilder)
    {
        this(configBuilder.getInteger(GERMLINE_PON_MARGIN), configBuilder.getValue(GERMLINE_PON_BED_SV_FILE),
                configBuilder.getValue(GERMLINE_PON_BED_SGL_FILE), false);

        buildSpecificSglRegion(RefGenomeVersion.from(configBuilder));
    }

    public PonCache(final int margin, final String ponSvFile, final String ponSglFile, boolean allowUnordered)
    {
        mSvRegions = Maps.newHashMap();
        mSglRegions = Maps.newHashMap();
        mAllowUnordered = allowUnordered;
        mHasValidData = true;

        mSpecificSglFusionRegions = Lists.newArrayList();

        mPositionMargin = margin;

        if(ponSvFile != null)
            loadPonSvFile(ponSvFile);

        if(ponSglFile != null)
            loadPonSglFile(ponSglFile);

        mCurrentSglIndex = 0;
        mCurrentSvIndex = 0;
    }

    private void buildSpecificSglRegion(final RefGenomeVersion refGenomeVersion)
    {
        if(refGenomeVersion.is37())
        {
            mSpecificSglFusionRegions.addAll(MULTI_MAPPED_ALT_ALIGNMENT_REGIONS_V37);
            mSpecificSglFusionRegions.addAll(SSX2_REGIONS_V37);
        }
        else
        {
            mSpecificSglFusionRegions.addAll(MULTI_MAPPED_ALT_ALIGNMENT_REGIONS_V38);
            mSpecificSglFusionRegions.addAll(SSX2_REGIONS_V38);
        }
    }

    public boolean hasValidData() { return mHasValidData; }

    public Map<String,List<PonSvRegion>> svRegions() { return mSvRegions; }
    public Map<String,List<PonSglRegion>> sglRegions() { return mSglRegions; }

    public void annotateVariants(final Map<String,List<Breakend>> chrBreakendMap)
    {
        for(Map.Entry<String,List<Breakend>> entry : chrBreakendMap.entrySet())
        {
            String chromosome = entry.getKey();

            mCurrentSvIndex = 0;
            mCurrentSglIndex = 0;

            List<PonSglRegion> sglRegions = mSglRegions.get(chromosome);
            List<PonSvRegion> svRegions = mSvRegions.get(chromosome);

            int lastPosStart = -1;

            for(Breakend breakend : entry.getValue())
            {
                if(!breakend.isSgl() && breakend.isEnd()) // only looked up on the first breakend
                    continue;

                if(lastPosStart > 0 && breakend.Position < lastPosStart)
                {
                    SV_LOGGER.error("var({}) out of order");
                }

                lastPosStart = breakend.Position;

                Variant var = breakend.sv();

                if(var.ponCount() > 0) // ignore if already annotated
                    continue;

                int ponCount = 0;
                int compareCount = 0;

                if(breakend.isSgl())
                {
                    ponCount = findSglPonMatch(sglRegions, var);

                    if(mRunChecks)
                        compareCount = findSglPonMatchBasic(sglRegions, var);
                }
                else
                {
                    ponCount = findSvPonMatch(svRegions, var);

                    if(mRunChecks)
                        compareCount = findSvPonMatchBasic(svRegions, var);
                }

                if(mRunChecks && ponCount != compareCount)
                {
                    SV_LOGGER.error("var({}) incorrect pon({}) vs basic({})", var, ponCount, compareCount);
                }

                if(ponCount > 0)
                {
                    var.setPonCount(ponCount);
                    var.addFilter(PON);
                }
            }
        }
    }

    private static int[] breakendMargin(final Breakend breakend)
    {
        int[] margins = new int[SE_PAIR];
        int inexactHomology = abs(breakend.IsStart ? breakend.InexactHomology.Start : breakend.InexactHomology.End);
        margins[SE_START] = min(breakend.ConfidenceInterval.Start, -inexactHomology);
        margins[SE_END] = max(breakend.ConfidenceInterval.End, inexactHomology);
        return margins;
    }

    private int findSvPonMatch(final List<PonSvRegion> regions, final Variant var)
    {
        if(regions == null)
            return 0;

        final int[] marginStart = breakendMargin(var.breakendStart());
        final int[] marginEnd = breakendMargin(var.breakendEnd());

        BaseRegion svStart = new BaseRegion(
                var.posStart() + marginStart[SE_START] - mPositionMargin,
                var.posStart() + marginStart[SE_END] + mPositionMargin);

        for(; mCurrentSvIndex < regions.size(); ++mCurrentSvIndex)
        {
            PonSvRegion region = regions.get(mCurrentSvIndex);

            if(region.overlapsStart(svStart))
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
        int maxPonCount = 0;
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

                if(region.matches(svStart, svEnd, var.orientStart(), var.orientEnd()) && region.PonCount > maxPonCount)
                    maxPonCount = region.PonCount;

                if(searchUp)
                    ++currentIndex;
                else
                    --currentIndex;
            }
        }

        return maxPonCount;
    }

    private boolean matchesSpecificSglFusionRegion(final Variant var)
    {
        // ignore if the SGL has an alt-mapping in a specific known fusion & poorly mapped region
        String alignmentsStr = var.breakendStart().Context.getAttributeAsString(INSALN, "");

        if(alignmentsStr.isEmpty())
            return false;

        List<AlternativeAlignment> alignments = AlternativeAlignment.fromVcfTag(alignmentsStr);

        if(alignments == null)
            return false;

        for(AlternativeAlignment altAlignment : alignments)
        {
            if(mSpecificSglFusionRegions.stream().anyMatch(x -> x.containsPosition(altAlignment.Chromosome, altAlignment.Position)))
                return true;
        }

        return false;
    }

    private int findSglPonMatch(final List<PonSglRegion> regions, final Variant var)
    {
        // ignore if the SGL has an alt-mapping in a specific known fusion & poorly mapped region
        if(matchesSpecificSglFusionRegion(var))
            return 0;

        if(regions == null)
            return 0;

        final int[] marginStart = breakendMargin(var.breakendStart());

        BaseRegion svStart = new BaseRegion(
                var.posStart() + marginStart[SE_START] - mPositionMargin,
                var.posStart() + marginStart[SE_END] + mPositionMargin);

        for(; mCurrentSglIndex < regions.size(); ++mCurrentSglIndex)
        {
            PonSglRegion region = regions.get(mCurrentSglIndex);

            if(region.overlaps(svStart))
            {
                // test the PON entries around this position
                return findSglPonMatch(regions, var, svStart, mCurrentSglIndex);
            }

            // exit if the PON is now past this point then retreat one position
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

        int maxPonCount = 0;
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

                if(region.matches(svStart, var.orientStart()) && region.PonCount > maxPonCount)
                    maxPonCount = region.PonCount;

                if(searchUp)
                    ++currentIndex;
                else
                    --currentIndex;
            }
        }

        return maxPonCount;
    }

    private int findSvPonMatchBasic(final List<PonSvRegion> regions, final Variant var)
    {
        final int[] marginStart = breakendMargin(var.breakendStart());
        final int[] marginEnd = breakendMargin(var.breakendEnd());

        BaseRegion svStart = new BaseRegion(
                var.posStart() + marginStart[SE_START] - mPositionMargin,
                var.posStart() + marginStart[SE_END] + mPositionMargin);

        for(int i = 0; i < regions.size(); ++i)
        {
            PonSvRegion region = regions.get(i);

            if(region.overlapsStart(svStart))
            {
                // test the PON entries around this position
                ChrBaseRegion svEnd = new ChrBaseRegion(
                        var.chromosomeEnd(),
                        var.posEnd() + marginEnd[SE_START] - mPositionMargin,
                        var.posEnd() + marginEnd[SE_END] + mPositionMargin);

                return findPonMatch(regions, var, svStart, svEnd, i);
            }

            if(region.RegionStart.start() > svStart.end())
                break;
        }

        return 0;
    }

    private int findSglPonMatchBasic(final List<PonSglRegion> regions, final Variant var)
    {
        final int[] marginStart = breakendMargin(var.breakendStart());

        BaseRegion svStart = new BaseRegion(
                var.posStart() + marginStart[SE_START] - mPositionMargin,
                var.posStart() + marginStart[SE_END] + mPositionMargin);

        for(int i = 0; i < regions.size(); ++i)
        {
            PonSglRegion region = regions.get(i);

            if(region.overlaps(svStart))
            {
                // test the PON entries around this position
                return findSglPonMatch(regions, var, svStart, 0);
            }

            // exit if the PON is now past this point then retreat one position
            if(region.Region.start() > svStart.end())
                break;
        }

        return 0;
    }

    public void checkSorted()
    {
        for(List<PonSvRegion> regions : mSvRegions.values())
        {
            Collections.sort(regions);
        }

        for(List<PonSglRegion> regions : mSglRegions.values())
        {
            Collections.sort(regions);
        }
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
            ChrBaseRegion lastRegionStart = null;

            while((line = fileReader.readLine()) != null)
            {
                PonSvRegion ponRegion = PonSvRegion.fromBedRecord(line);

                if(!ponRegion.RegionStart.Chromosome.equals(currentChr))
                {
                    currentChr = ponRegion.RegionStart.Chromosome;

                    if(!mSvRegions.containsKey(ponRegion.RegionStart.Chromosome))
                    {
                        svRegions = Lists.newArrayList();
                        mSvRegions.put(ponRegion.RegionStart.Chromosome, svRegions);
                    }
                    else
                    {
                        svRegions = mSvRegions.get(ponRegion.RegionStart.Chromosome);
                    }

                    lastRegionStart = null;
                }

                svRegions.add(ponRegion);
                ++itemCount;

                if(!mAllowUnordered && lastRegionStart != null && lastRegionStart.start() > ponRegion.RegionStart.start())
                {
                    SV_LOGGER.warn("SV PON not ordered: last({}) vs this({})", lastRegionStart, ponRegion.RegionStart);
                    mHasValidData = false;
                }

                lastRegionStart = ponRegion.RegionStart;
            }

            SV_LOGGER.info("loaded {} SV PON records from file({})", itemCount, filename);
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to load SV PON file({}): {}", filename, e.toString());
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
            ChrBaseRegion lastRegion = null;

            while((line = fileReader.readLine()) != null)
            {
                PonSglRegion ponRegion = PonSglRegion.fromBedRecord(line);

                if(!ponRegion.Region.Chromosome.equals(currentChr))
                {
                    currentChr = ponRegion.Region.Chromosome;

                    if(!mSglRegions.containsKey(ponRegion.Region.Chromosome))
                    {
                        sglRegions = Lists.newArrayList();
                        mSglRegions.put(ponRegion.Region.Chromosome, sglRegions);
                    }
                    else
                    {
                        sglRegions = mSglRegions.get(ponRegion.Region.Chromosome);
                    }

                    lastRegion = null;
                }

                sglRegions.add(ponRegion);
                ++itemCount;
                
                if(!mAllowUnordered && lastRegion != null && lastRegion.start() > ponRegion.Region.start())
                {
                    SV_LOGGER.warn("SGL PON not ordered: last({}) vs this({})", lastRegion, ponRegion.Region);
                    mHasValidData = false;
                }

                lastRegion = ponRegion.Region;
            }

            SV_LOGGER.info("loaded {} SGL PON records from file({})", itemCount, filename);
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to load SGL PON file({}): {}", filename, e.toString());
            mHasValidData = false;
        }
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

    /*
    @VisibleForTesting
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

    @VisibleForTesting
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
    */
}
