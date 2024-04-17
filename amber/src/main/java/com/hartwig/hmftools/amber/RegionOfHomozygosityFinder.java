package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates.COORDS_37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates.COORDS_38;
import static com.hartwig.hmftools.common.utils.collection.Multimaps.filterEntries;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/**
 * Regions of Homozygosity(ROH) Algorithm
 *
 * Amber outputs a file which contains continuous regions of non high confidence heterozygous points.
 * A region must meet the following criteria to be considered a region of homozygosity:
 * Region must have at least 50 consecutive SNP locations
 * No 5 of any consecutive 50 locations heterozygous
 * Region must be at least 500,000 bases in length
 *
 **/
public class RegionOfHomozygosityFinder
{
    enum Zygosity
    {
        HETEROZYGOUS,
        HOMOZYGOUS,
        UNCLEAR
    }

    static class LocusZygosity
    {
        public final int position;
        public final Zygosity zygosity;

        LocusZygosity(int pos, Zygosity zygosity)
        {
            this.position = pos;
            this.zygosity = zygosity;
        }
    }

    private static final int EXCLUDED_REGION_EXPAND_SIZE = 1_000_000;
    private static final int HALF_CENTROMERE_SIZE = 1_500_000;

    private final List<GenomeRegion> mExcludedRegions;
    private final RefGenomeVersion mRefGenomeVersion;
    private final double mMinDepthPercent;
    private final double mMaxDepthPercent;
    private final int mMinHomozygousRegionSize;
    private final int mMinSnpLociCount;
    private final int mSiteWindowSize;
    private final int mMaxHetInWindow;

    public RegionOfHomozygosityFinder(
            final RefGenomeVersion refGenomeVersion, double minDepthPercent, double maxDepthPercent,
            int minHomozygousRegionSize, int minSnpLociCount, int siteWindowSize, int maxHetInWindow) throws IOException
    {
        mRefGenomeVersion = refGenomeVersion;
        mMinDepthPercent = minDepthPercent;
        mMaxDepthPercent = maxDepthPercent;
        mMinHomozygousRegionSize = minHomozygousRegionSize;
        mMinSnpLociCount = minSnpLociCount;
        mSiteWindowSize = siteWindowSize;
        mMaxHetInWindow = maxHetInWindow;
        mExcludedRegions = Lists.newArrayList();

        loadExcludedRegions();
    }

    public RegionOfHomozygosityFinder(RefGenomeVersion refGenomeVersion, double minDepthPercent, double maxDepthPercent) throws IOException
    {
        this(refGenomeVersion,
                minDepthPercent, maxDepthPercent,
                AmberConstants.HOMOZYGOUS_REGION_MIN_SIZE,
                AmberConstants.HOMOZYGOUS_REGION_MIN_SNP_LOCI_COUNT,
                AmberConstants.HOMOZYGOUS_REGION_WINDOW_SIZE, AmberConstants.HOMOZYGOUS_REGION_MAX_HET_IN_WINDOW);
    }

    public List<RegionOfHomozygosity> findRegions(final ListMultimap<Chromosome, PositionEvidence> baseDepths)
    {
        final Predicate<PositionEvidence> depthFilter = new BaseDepthFilter(mMinDepthPercent, mMaxDepthPercent, baseDepths);
        final ListMultimap<Chromosome, PositionEvidence> filteredBaseDepths = filterEntries(baseDepths, depthFilter);

        List<RegionOfHomozygosity> homozygousRegions = new ArrayList<RegionOfHomozygosity>();

        Set<Chromosome> chromosomeSet = baseDepths.keySet();

        for(Chromosome chromosome : chromosomeSet)
        {
            if(chromosome.isAllosome())
                continue;
            
            homozygousRegions.addAll(findRegionsForChromosome(chromosome, toLocusZygosityList(filteredBaseDepths.get(chromosome))));
        }

        Collections.sort(homozygousRegions);

        return homozygousRegions;
    }

    public List<RegionOfHomozygosity> findRegionsForChromosome(Chromosome chromosome, final List<LocusZygosity> bafSites)
    {
        var homozygousRegions = new ArrayList<RegionOfHomozygosity>();

        // keep track of where the sliding window is
        // following are indices
        int slidingWindowStartIndex = -1;
        int rohStartIndex = -1;
        int rohEndIndex = -1;

        int numHetInWindow = 0;

        for(int i = 0; i < bafSites.size(); ++i)
        {
            LocusZygosity bafSite = bafSites.get(i);
            int position = bafSite.position;
            Zygosity zygosity = bafSite.zygosity;

            if(rohStartIndex != -1)
            {
                boolean crossedExcludedRegion = overlapWithExcludedRegions(chromosome.toString(),
                        rohEndIndex != -1 ? bafSites.get(rohEndIndex).position : position, position);

                // first thing we need to check if the previous region is ending
                if(zygosity == Zygosity.HETEROZYGOUS)
                    ++numHetInWindow;

                // say sliding window size is 3, and i = 4. we want to remove i = 4 - 3 = 1
                assert (slidingWindowStartIndex >= i - mSiteWindowSize);

                if(slidingWindowStartIndex <= i - mSiteWindowSize)
                {
                    // one location has gone out of the sliding window
                    if(bafSites.get(slidingWindowStartIndex).zygosity == Zygosity.HETEROZYGOUS)
                    {
                        // if it is het, then we reduce the count
                        --numHetInWindow;
                    }
                    ++slidingWindowStartIndex;
                }

                if(numHetInWindow > mMaxHetInWindow || crossedExcludedRegion)
                {
                    // this stretch of homozygosity ends
                    RegionOfHomozygosity region = createRegionIfPassFilter(bafSites, chromosome, rohStartIndex, rohEndIndex);
                    if(region != null)
                    {
                        homozygousRegions.add(region);
                    }
                    rohStartIndex = -1;
                    rohEndIndex = -1;
                }
            }

            if(zygosity == Zygosity.HOMOZYGOUS && !overlapWithExcludedRegions(chromosome.toString(), position, position))
            {
                if(rohStartIndex != -1)
                {
                    // we are continuing a region, so update the end
                    rohEndIndex = i;
                }
                else
                {
                    // we can start a new region here
                    rohStartIndex = i;
                    rohEndIndex = i;

                    // we start the window again
                    slidingWindowStartIndex = i;
                    numHetInWindow = 0;
                }
            }
        }

        // this stretch of homozygosity ends
        RegionOfHomozygosity region = createRegionIfPassFilter(bafSites, chromosome, rohStartIndex, rohEndIndex);
        if(region != null)
        {
            // this is a region
            homozygousRegions.add(region);
        }

        return homozygousRegions;
    }

    // filter them by following rules:
    // 1. Do not allow regions to pass over the centromere or the adjacent heterochromatin segments on chr 1 and 9 (chr1:124535434-142535434;
    //    chr16:38335801-46335801, chr9:50367679-65367679).  Also ignore any BAF points within 1M bases of these gaps
    // 2. Do not analyse chromsome X if amber determines the gender as male.
    // 3. Exclude the p arm of chromosome 13,14,15,21 and 22.
    //

    public void loadExcludedRegions() throws IOException
    {
        String resourcePath = null;
        switch(mRefGenomeVersion)
        {
            case V37:
                resourcePath = "rohExcluded.37.bed";
                break;
            case V38:
                resourcePath = "rohExcluded.38.bed";
                break;
        }

        List<ChrBaseRegion> excludedRegions = AmberUtils.loadBedFromResource(resourcePath);

        // add a buffer around each centromere
        RefGenomeCoordinates refCcoord = mRefGenomeVersion == RefGenomeVersion.V37 ? COORDS_37 : COORDS_38;

        for(Map.Entry<Chromosome, Integer> entry : refCcoord.Centromeres.entrySet())
        {
            int centromerePosition = entry.getValue();

            excludedRegions.add(new ChrBaseRegion(
                    entry.getKey().toString(),
                    centromerePosition - HALF_CENTROMERE_SIZE,
                    centromerePosition + HALF_CENTROMERE_SIZE));
        }

        excludedRegions.forEach(x -> mExcludedRegions.add(x.genomeRegion()));
    }

    private boolean overlapWithExcludedRegions(String chromosome, int position1, int position2)
    {
        assert(position1 <= position2);
        for(GenomeRegion genomeRegion: mExcludedRegions)
        {
            // we want to check if we overlap
            if(genomeRegion.chromosome().equals(chromosome) &&
                    position1 <= genomeRegion.end() + EXCLUDED_REGION_EXPAND_SIZE &&
                    position2 >= genomeRegion.start() - EXCLUDED_REGION_EXPAND_SIZE)
            {
                return true;
            }
        }
        return false;
    }

    // we apply some filters
    @Nullable
    private RegionOfHomozygosity createRegionIfPassFilter(
            @NotNull List<LocusZygosity> bafSites,
            Chromosome chromosome, int regionStartIndex, int regionEndIndex)
    {
        if(regionStartIndex <= 0 || regionEndIndex <= 0)
            return null;

        int numHomozygous = 0;
        int numHeterozygous = 0;
        int numUnclear = 0;

        // we want to count
        for(int i = regionStartIndex; i <= regionEndIndex; ++i)
        {
            Zygosity zygosity = bafSites.get(i).zygosity;
            switch (zygosity)
            {
                case HOMOZYGOUS:
                    ++numHomozygous;
                    break;
                case HETEROZYGOUS:
                    ++numHeterozygous;
                    break;
            }
        }

        int snpCount = numHomozygous + numHeterozygous + numUnclear;
        if(snpCount < mMinSnpLociCount)
        {
            return null;
        }
        int regionStart = bafSites.get(regionStartIndex).position;
        int regionEnd = bafSites.get(regionEndIndex).position;

        if((regionEnd - regionStart) >= mMinHomozygousRegionSize)
        {
            // this is a region
            return new RegionOfHomozygosity(chromosome, regionStart, regionEnd, numHomozygous, numHeterozygous, numUnclear);
        }
        return null;
    }

    // we use same definition as purple here.
    static boolean isAlleleHomozygous(int totalCount, int alleleCount) {
        if(totalCount == alleleCount)
        {
            return true;
        }

        boolean isHighVaf = alleleCount > 0.75 * totalCount;
        if(isHighVaf)
        {
            double p = new PoissonDistribution(totalCount / 2d).cumulativeProbability(totalCount - alleleCount);
            return p < 0.005;
        }

        return false;
    }

    static Zygosity calcZygosity(PositionEvidence baseDepth)
    {
        if(isAlleleHomozygous(baseDepth.ReadDepth, baseDepth.RefSupport) || isAlleleHomozygous(baseDepth.ReadDepth, baseDepth.AltSupport))
        {
            return Zygosity.HOMOZYGOUS;
        }
        return Zygosity.HETEROZYGOUS;
    }

    static LocusZygosity toLocusZygosity(PositionEvidence depth)
    {
        return new LocusZygosity(depth.position(), calcZygosity(depth));
    }

    private static List<LocusZygosity> toLocusZygosityList(@NotNull final List<PositionEvidence> bafs)
    {
        var locusZygosityList = new ArrayList<LocusZygosity>();

        for(var depth : bafs)
        {
            locusZygosityList.add(toLocusZygosity(depth));
        }

        return locusZygosityList;
    }
}
