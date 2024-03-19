package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.metrics.BamMetricsSummary.BAM_METRICS_FILE_ID;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;

public class OffTargetFragments
{
    private final List<Fragment> mFragments;
    private final List<EnrichedSite> mEnrichedSites;

    private int mCurrentLowerPosition;
    private int mCurrentUpperPosition;

    private final Map<Integer,Integer> mOverlapCounts;

    private final int mFragmentOverlapWriteThreshold;

    public OffTargetFragments(final int fragmentOverlapWriteThreshold)
    {
        mFragmentOverlapWriteThreshold = fragmentOverlapWriteThreshold;
        mFragments = Lists.newArrayList();
        mEnrichedSites = Lists.newArrayList();
        mCurrentLowerPosition = 0;
        mCurrentUpperPosition = 0;
        mOverlapCounts = Maps.newHashMap();
    }

    public Map<Integer,Integer> fragmentOverlapCounts() { return mOverlapCounts; }

    public void addRead(final SAMRecord read, final boolean isLocalConcordantFragment, final int readMateEnd)
    {
        if(read.getReadUnmappedFlag())
            return;

        if(mFragments.stream().anyMatch(x -> x.ReadId.equals(read.getReadName())))
            return;

        Fragment newFragment = new Fragment(read, isLocalConcordantFragment, readMateEnd);

        if(newFragment.AlignmentStart > mCurrentUpperPosition)
            processOverlappingFragments();

        if(mFragments.isEmpty())
            mCurrentLowerPosition = newFragment.AlignmentStart;

        mCurrentUpperPosition = max(mCurrentUpperPosition, newFragment.AlignmentEnd);
        mFragments.add(newFragment);
    }

    private void processOverlappingFragments()
    {
        int maxOverlaps = 0;
        Fragment peakOverlapFragment = null;

        for(Fragment first : mFragments)
        {
            int overlaps = 1; // count itself since this value is a measure of depth at this location

            for(Fragment second : mFragments)
            {
                if(first == second)
                    continue;

                if(positionsOverlap(first.AlignmentStart, first.AlignmentEnd, second.AlignmentStart, second.AlignmentEnd))
                    ++overlaps;
                else if(second.AlignmentStart > first.AlignmentEnd) // all subsequent fragments will also not overlap
                    break;
            }

            Integer frequency = mOverlapCounts.get(overlaps);
            mOverlapCounts.put(overlaps, frequency != null ? frequency + 1 : 1);

            if(overlaps > maxOverlaps)
            {
                maxOverlaps = overlaps;
                peakOverlapFragment = first;
            }
        }

        if(mFragmentOverlapWriteThreshold > 0 && maxOverlaps >= mFragmentOverlapWriteThreshold)
        {
            mEnrichedSites.add(new EnrichedSite(
                    new ChrBaseRegion(mFragments.get(0).Chromosome, mCurrentLowerPosition, mCurrentUpperPosition),
                    new BaseRegion(peakOverlapFragment.AlignmentStart, peakOverlapFragment.AlignmentEnd),
                    mFragments.size(), maxOverlaps));
        }

        mFragments.clear();
    }

    private class Fragment
    {
        public final String Chromosome;
        public final String ReadId;
        public final int AlignmentStart;
        public final int AlignmentEnd;
        public final boolean Local;

        public Fragment(final SAMRecord read, final boolean isLocalConcordantFragment, final int readMateEnd)
        {
            Chromosome = read.getReferenceName();
            ReadId = read.getReadName();
            Local = isLocalConcordantFragment;

            if(!isLocalConcordantFragment)
            {
                AlignmentStart = read.getAlignmentStart();
                AlignmentEnd = read.getAlignmentEnd();
            }
            else
            {
                AlignmentStart = min(read.getAlignmentStart(), read.getMateAlignmentStart());
                AlignmentEnd = max(read.getAlignmentEnd(), readMateEnd);
            }
        }

        public String toString()
        {
            return format("coords(%s:%d-%d) %s id(%s)",
                    Chromosome, AlignmentStart, AlignmentEnd, Local ? "local" : "discordant", ReadId);
        }
    }

    private class EnrichedSite
    {
        public final ChrBaseRegion Region;
        public final BaseRegion PeakRegion;
        public final int Fragments;
        public final int MaxOverlap;

        public EnrichedSite(final ChrBaseRegion region, final BaseRegion peakRegion, final int fragments, final int maxOverlap)
        {
            Region = region;
            PeakRegion = peakRegion;
            Fragments = fragments;
            MaxOverlap = maxOverlap;
        }

        public String toString()
        {
            return format("region(%s) fragments(%d) maxOverlaps(%d)", Region, Fragments, MaxOverlap);
        }
    }

    public List<String> enrichedFragmentSites()
    {
        if(mEnrichedSites.isEmpty())
            return Collections.emptyList();

        List<String> enrichedDataInfo = Lists.newArrayListWithCapacity(mEnrichedSites.size());

        for(EnrichedSite enrichedSite : mEnrichedSites)
        {
            enrichedDataInfo.add(format("%s\t%s\t%d\t%d\t%d\t%d\t%d",
                    enrichedSite.Region.Chromosome, enrichedSite.Region.start(), enrichedSite.Region.end(),
                    enrichedSite.PeakRegion.start(), enrichedSite.PeakRegion.end(),
                    enrichedSite.Fragments, enrichedSite.MaxOverlap));
        }

        mEnrichedSites.clear();

        return enrichedDataInfo;
    }

    public String toString() { return format("fragments(%d) range(%d-%d)", mFragments.size(), mCurrentLowerPosition, mCurrentUpperPosition); }

    public static void writeOverlapCounts(final MetricsConfig config, final Map<Integer,Integer> offTargetOverlapCounts)
    {
        try
        {
            String filename = config.formFilename("off_target_overlaps");

            BufferedWriter writer = createBufferedWriter(filename, false);
            writer.write("FragmentOverlaps\tFrequency");
            writer.newLine();

            List<Integer> overlaps = offTargetOverlapCounts.keySet().stream().collect(Collectors.toList());
            Collections.sort(overlaps);

            for(Integer overlapCount : overlaps)
            {
                writer.write(format("%d\t%d", overlapCount, offTargetOverlapCounts.get(overlapCount)));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to writer off-target fragment overlaps: {}", e.toString());
        }
    }

    public static BufferedWriter initialiseEnrichedRegionWriter(final MetricsConfig config)
    {
        try
        {
            String filename = config.formFilename("off_target_enriched");

            BufferedWriter writer = createBufferedWriter(filename, false);
            writer.write("Chromosome\tPosStart\tPosEnd\tPeakPosStart\tPeakPosEnd\tFragments\tMaxOverlap");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to initialise high fragment overlaps writer: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeEnrichedRegions(final BufferedWriter writer, final List<String> fragmentData)
    {
        try
        {
            for(String fragment : fragmentData)
            {
                writer.write(fragment);
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write high fragment overlaps writer: {}", e.toString());
        }
    }
}
