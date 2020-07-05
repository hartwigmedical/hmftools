package com.hartwig.hmftools.svtools.cohort;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionWithin;
import static com.hartwig.hmftools.svtools.cohort.LineElementType.NONE;
import static com.hartwig.hmftools.svtools.cohort.LineElementType.moreKnown;

import java.util.List;

import com.hartwig.hmftools.common.utils.sv.SvRegion;

import org.apache.commons.compress.utils.Lists;

public class LineClusterData
{
    public final String SampleId;

    public final int ClusterId;

    public final List<LineRegion> LineRegions;

    public final List<SvRegion> InsertedRegions;

    public final List<LineClusterData> MatchedClusters;

    public LineClusterData(final String sampleId, final int clusterId, final SvRegion lineRegion, final LineElementType lineType)
    {
        SampleId = sampleId;
        ClusterId = clusterId;

        LineRegions = Lists.newArrayList();
        LineRegions.add(new LineRegion(lineRegion, lineType));

        InsertedRegions = Lists.newArrayList();
        MatchedClusters = Lists.newArrayList();
    }

    public static int PROXIMATE_DISTANCE = 10000;

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

    public static boolean regionWithinBounds(final SvRegion region, final String chromosome, final int position)
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

        SvRegion newRegion = new SvRegion(chromosome, position, position);
        LineRegions.add(new LineRegion(newRegion, lineType));
    }

    public void addInsertData(final String chromosome, final int position)
    {
        if(chromosome.equals("0"))
            return;

        for(SvRegion insertRegion : InsertedRegions)
        {
            if(regionWithinBounds(insertRegion, chromosome, position))
            {
                insertRegion.setStart(min(insertRegion.start(), position));
                insertRegion.setEnd(max(insertRegion.end(), position));
                return;
            }
        }

        InsertedRegions.add(new SvRegion(chromosome, position, position));
    }

    public boolean matches(final LineClusterData other)
    {
        return LineRegions.stream().anyMatch(x -> other.LineRegions.stream().anyMatch(y -> y.Region.overlaps(x.Region)));
    }

    public int sampleCount() { return 1 + MatchedClusters.size(); }
    public int insertRegionsCount()
    {
        int matchedInsertTotal = MatchedClusters.stream().mapToInt(x -> x.InsertedRegions.size()).sum();
        return InsertedRegions.size() + matchedInsertTotal; }

    public LineRegion primaryRegion()
    {
        LineRegion maxRegion = LineRegions.get(0);

        for(int i = 1; i < LineRegions.size(); ++i)
        {
            if(LineRegions.get(i).BreakendCount > maxRegion.BreakendCount)
            {
                maxRegion = LineRegions.get(i);
            }
        }

        return maxRegion;
    }

    public String toString()
    {
        final LineRegion primaryRegion = primaryRegion();
        return String.format("source(%s) type(%s) inserts(%d) sourceLocations(%d) sample(%s:%s)",
                primaryRegion, primaryRegion.LineType, InsertedRegions.size(), LineRegions.size(),
            SampleId, ClusterId);
    }

    public String sampleClusterData()
    {
        return String.format("%s,%d,%d", SampleId, ClusterId, InsertedRegions.size());
    }

}
