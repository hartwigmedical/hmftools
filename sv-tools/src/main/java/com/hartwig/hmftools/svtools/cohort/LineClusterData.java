package com.hartwig.hmftools.svtools.cohort;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.svtools.cohort.LineElementType.KNOWN;
import static com.hartwig.hmftools.svtools.cohort.LineElementType.KNOWN_SUSPECT;
import static com.hartwig.hmftools.svtools.cohort.LineElementType.NONE;
import static com.hartwig.hmftools.svtools.cohort.LineElementType.SUSPECT;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.compress.utils.Lists;

public class LineClusterData
{
    public final String SampleId;

    public final int ClusterId;

    public final List<LineRegion> LineRegions;

    public final List<ChrBaseRegion> InsertedRegions;

    public final List<LineClusterData> MatchedClusters;

    private LineRegion mPrimaryRegion;

    public LineClusterData(final String sampleId, final int clusterId, final ChrBaseRegion lineRegion, final LineElementType lineType)
    {
        SampleId = sampleId;
        ClusterId = clusterId;

        LineRegions = Lists.newArrayList();
        LineRegions.add(new LineRegion(lineRegion, lineType));
        mPrimaryRegion = LineRegions.get(0);

        InsertedRegions = Lists.newArrayList();
        MatchedClusters = Lists.newArrayList();
    }

    public static final int PROXIMATE_DISTANCE = 10000;

    public void addSvData(final String[] chromosomes, final int[] positions, final LineElementType[] lineTypes)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(chromosomes.equals("0"))
                continue;

            if(lineTypes[se] == NONE)
            {
                addInsertData(chromosomes[se], positions[se]);
            }
            else
            {
                addSourceData(chromosomes[se], positions[se], lineTypes[se]);
            }
        }
    }

    private static boolean regionWithinBounds(final ChrBaseRegion region, final String chromosome, final int position)
    {
        return region.Chromosome.equals(chromosome) && positionWithin(
                position, region.start() - PROXIMATE_DISTANCE, region.end() + PROXIMATE_DISTANCE);
    }

    public void addSourceData(final String chromosome, final int position, final LineElementType lineType)
    {
        if(chromosome.equals("0"))
            return;

        LineRegion matchedRegion = LineRegions.stream()
                .filter(x -> regionWithinBounds(x.Region, chromosome, position)).findFirst().orElse(null);

        if(matchedRegion != null)
        {
            ++matchedRegion.BreakendCount;
            matchedRegion.Region.setStart(min(matchedRegion.Region.start(), position));
            matchedRegion.Region.setEnd(max(matchedRegion.Region.end(), position));
            return;
        }

        ChrBaseRegion newRegion = new ChrBaseRegion(chromosome, position, position);
        LineRegions.add(new LineRegion(newRegion, lineType));
        setPrimaryRegion();
    }

    public void addInsertData(final String chromosome, final int position)
    {
        if(chromosome.equals("0"))
            return;

        for(ChrBaseRegion insertRegion : InsertedRegions)
        {
            if(regionWithinBounds(insertRegion, chromosome, position))
            {
                insertRegion.setStart(min(insertRegion.start(), position));
                insertRegion.setEnd(max(insertRegion.end(), position));
                return;
            }
        }

        InsertedRegions.add(new ChrBaseRegion(chromosome, position, position));
    }

    public boolean matches(final LineClusterData other)
    {
        if(!mPrimaryRegion.Region.Chromosome.equals(other.primaryRegion().Region.Chromosome))
            return false;

        return positionsOverlap(
                mPrimaryRegion.Region.start() - PROXIMATE_DISTANCE, mPrimaryRegion.Region.end() + PROXIMATE_DISTANCE,
        other.primaryRegion().Region.start(), other.primaryRegion().Region.end());
    }

    public int sampleCount() { return 1 + MatchedClusters.size(); }
    public int insertRegionsCount()
    {
        int matchedInsertTotal = MatchedClusters.stream().mapToInt(x -> x.InsertedRegions.size()).sum();
        return InsertedRegions.size() + matchedInsertTotal;
    }

    public ChrBaseRegion getCombinedPrimarySourceRegion()
    {
        if(MatchedClusters.isEmpty())
            return mPrimaryRegion.Region;

        ChrBaseRegion combinedRegion = new ChrBaseRegion(mPrimaryRegion.Region.Chromosome, mPrimaryRegion.Region.Positions);

        for(LineClusterData otherCluster : MatchedClusters)
        {
            combinedRegion.setStart(min(combinedRegion.start(), otherCluster.primaryRegion().Region.start()));
            combinedRegion.setEnd(max(combinedRegion.end(), otherCluster.primaryRegion().Region.end()));
        }

        return combinedRegion;
    }

    public LineRegion primaryRegion() { return mPrimaryRegion; }

    private void setPrimaryRegion()
    {
        // take the highest of line type followed by highest breakend count
        List<LineRegion> lineRegions = LineRegions.stream().filter(x -> x.LineType == KNOWN).collect(Collectors.toList());

        if(lineRegions.isEmpty())
            lineRegions = LineRegions.stream().filter(x -> x.LineType == KNOWN_SUSPECT).collect(Collectors.toList());

        if(lineRegions.isEmpty())
            lineRegions = LineRegions.stream().filter(x -> x.LineType == SUSPECT).collect(Collectors.toList());

        LineRegion maxRegion = lineRegions.get(0);

        for(int i = 1; i < lineRegions.size(); ++i)
        {
            if(lineRegions.get(i).BreakendCount > maxRegion.BreakendCount)
            {
                maxRegion = lineRegions.get(i);
            }
        }

        mPrimaryRegion = maxRegion;
    }

    public String toString()
    {
        return String.format("source(%s) type(%s) inserts(%d) sourceLocations(%d) sample(%s:%s)",
                mPrimaryRegion, mPrimaryRegion.LineType, InsertedRegions.size(), LineRegions.size(),
            SampleId, ClusterId);
    }

    public String sampleClusterData()
    {
        return String.format("%s,%d,%d,%d,%d,%d,%d",
                SampleId, ClusterId, mPrimaryRegion.Region.start(), mPrimaryRegion.Region.end(),
                LineRegions.size(), mPrimaryRegion.BreakendCount, InsertedRegions.size());
    }

}
