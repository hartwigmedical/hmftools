package com.hartwig.hmftools.linx.cohort;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHR_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHR_START;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENT_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENT_START;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;

public class PonBuilder
{
    // keyed by chrStart_orientStart_chrEnd_orientEnd,then posStart_posEnd for SVs
    private final Map<String,Map<String,Integer>> mSvPositionCountsMap;

    // keyed by chrStart_orientStart, then posStart for SGLs
    private final Map<String,Map<Integer,Integer>> mSglPositionCountsMap;

    public PonBuilder()
    {
        mSvPositionCountsMap = Maps.newHashMap();
        mSglPositionCountsMap = Maps.newHashMap();
    }

    public void registerSv(
            final String chrStart, final int posStart, final Orientation orientStart,
            final String chrEnd, final int posEnd, final Orientation orientEnd)
    {
        String chrOrientKey = formChrOrientKey(chrStart, orientStart, chrEnd, orientEnd);
        Map<String,Integer> positionCounts = mSvPositionCountsMap.get(chrOrientKey);

        if(positionCounts == null)
        {
            positionCounts = Maps.newHashMap();
            mSvPositionCountsMap.put(chrOrientKey, positionCounts);
        }

        String positionPair = formPositionPair(posStart, posEnd);
        Integer count = positionCounts.get(positionPair);
        positionCounts.put(positionPair, count != null ? count + 1 : 1);
    }

    public void registerSgl(final String chrStart, final int posStart, final Orientation orientStart)
    {
        String chrOrientKey = formChrOrientKey(chrStart, orientStart, null, null);
        Map<Integer,Integer> positionCounts = mSglPositionCountsMap.get(chrOrientKey);

        if(positionCounts == null)
        {
            positionCounts = Maps.newHashMap();
            mSglPositionCountsMap.put(chrOrientKey, positionCounts);
        }

        Integer count = positionCounts.get(posStart);
        positionCounts.put(posStart, count != null ? count + 1 : 1);
    }

    protected static String PAIR_DELIM = "_";

    protected static String formChrOrientKey(
            final String chrStart, final Orientation orientStart, @Nullable final String chrEnd, @Nullable final Orientation orientEnd)
    {
        if(chrEnd == null)
            return format("%s_%d", chrStart, orientStart.asByte());
        else
            return format("%s_%d_%s_%d", chrStart, orientStart.asByte(), chrEnd, orientEnd.asByte());
    }

    protected static String formPositionPair(final int posStart, final Integer posEnd)
    {
        return format("%d_%d", posStart, posEnd);
    }

    protected static final String FLD_COUNT = "Count";

    public void writePonCache(final CohortConfig config)
    {
        int ponMargin = config.PonMargin;

        try
        {
            String filename = config.formOutputFilename("pon_cache");
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            // definitional fields
            sj.add(FLD_CHR_START).add(FLD_POS_START).add(FLD_ORIENT_START);
            sj.add(FLD_CHR_END).add(FLD_POS_END).add(FLD_ORIENT_END);
            sj.add(FLD_COUNT);

            writer.write(sj.toString());

            writer.newLine();

            for(Map.Entry<String,Map<String,Integer>> chrEntry : mSvPositionCountsMap.entrySet())
            {
                String chrOrientKey = chrEntry.getKey();
                Map<String,Integer> locationMap = chrEntry.getValue();

                String[] parts = chrOrientKey.split(PAIR_DELIM, 4);
                String chrStart = parts[0];
                String orientStart = parts[1];
                String chrEnd = parts[2];
                String orientEnd = parts[3];

                List<PositionCount> positionPairs = Lists.newArrayList();

                for(Map.Entry<String,Integer> locationEntry : locationMap.entrySet())
                {
                    String positionKey = locationEntry.getKey();

                    String[] positionParts = positionKey.split(PAIR_DELIM, 2);

                    positionPairs.add(new PositionCount(
                            Integer.parseInt(positionParts[0]), Integer.parseInt(positionParts[1]), locationEntry.getValue()));
                }

                // merge proximate PON entries
                Collections.sort(positionPairs);

                int index = 0;
                while(index < positionPairs.size() - 1)
                {
                    PositionCount positionCount = positionPairs.get(index);

                    int nextIndex = index + 1;
                    while(nextIndex < positionPairs.size())
                    {
                        PositionCount nextPositionCount = positionPairs.get(nextIndex);

                        if(nextPositionCount.PosStartUpper - positionCount.PosStartLower <= ponMargin
                        && abs(nextPositionCount.PosEndUpper - positionCount.PosEndLower) <= ponMargin)
                        {
                            positionCount.PosStartUpper = max(positionCount.PosStartUpper, nextPositionCount.PosStartUpper);
                            int maxEnd = max(positionCount.PosEndUpper, nextPositionCount.PosEndUpper);
                            int minEnd = min(positionCount.PosEndUpper, nextPositionCount.PosEndUpper);
                            positionCount.PosEndLower = minEnd;
                            positionCount.PosEndUpper = maxEnd;
                            positionCount.Count += nextPositionCount.Count;
                            positionPairs.remove(nextIndex);
                        }
                        else
                        {
                            break;
                        }
                    }

                    ++index;
                }

                for(PositionCount positionCount : positionPairs)
                {
                    if(positionCount.Count <= 1)
                        continue;

                    sj = new StringJoiner(TSV_DELIM);

                    int posStart = (int)round((positionCount.PosStartLower + positionCount.PosStartUpper) * 0.5);
                    int posEnd = (int)round((positionCount.PosEndLower + positionCount.PosEndUpper) * 0.5);
                    sj.add(chrStart).add(String.valueOf(posStart)).add(orientStart);
                    sj.add(chrEnd).add(String.valueOf(posEnd)).add(orientEnd);

                    sj.add(String.valueOf(positionCount.Count));

                    writer.write(sj.toString());
                    writer.newLine();
                }
            }

            for(Map.Entry<String,Map<Integer,Integer>> chrEntry : mSglPositionCountsMap.entrySet())
            {
                String chrOrientKey = chrEntry.getKey();
                Map<Integer,Integer> locationMap = chrEntry.getValue();

                String[] parts = chrOrientKey.split(PAIR_DELIM, 2);
                String chrStart = parts[0];
                String orientStart = parts[1];

                List<PositionCount> positionPairs = Lists.newArrayList();

                for(Map.Entry<Integer,Integer> locationEntry : locationMap.entrySet())
                {
                    int position = locationEntry.getKey();
                    positionPairs.add(new PositionCount(position, 0, locationEntry.getValue()));
                }

                Collections.sort(positionPairs);

                int index = 0;
                while(index < positionPairs.size() - 1)
                {
                    PositionCount positionCount = positionPairs.get(index);

                    int nextIndex = index + 1;
                    while(nextIndex < positionPairs.size())
                    {
                        PositionCount nextPositionCount = positionPairs.get(nextIndex);

                        if(nextPositionCount.PosStartUpper - positionCount.PosStartLower <= ponMargin)
                        {
                            positionCount.PosStartUpper = max(positionCount.PosStartUpper, nextPositionCount.PosStartUpper);
                            positionCount.Count += nextPositionCount.Count;
                            positionPairs.remove(nextIndex);
                        }
                        else
                        {
                            break;
                        }
                    }

                    ++index;
                }

                for(PositionCount positionCount : positionPairs)
                {
                    if(positionCount.Count <= 1)
                        continue;

                    sj = new StringJoiner(TSV_DELIM);

                    int posStart = (int)round((positionCount.PosStartLower + positionCount.PosStartUpper) * 0.5);

                    sj.add(chrStart).add(String.valueOf(posStart)).add(orientStart);
                    sj.add("").add("").add("");

                    sj.add(String.valueOf(positionCount.Count));

                    writer.write(sj.toString());
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to write ponc cache: {}", e.toString());
        }
    }

    private class PositionCount implements Comparable<PositionCount>
    {
        public int PosStartLower;
        public int PosStartUpper;

        public int PosEndLower;
        public int PosEndUpper;

        public int Count;

        public PositionCount(final int posStart, final int posEnd, final int count)
        {
            PosStartLower = PosStartUpper = posStart;
            PosEndLower = PosEndUpper = posEnd;
            Count = count;
        }

        @Override
        public int compareTo(final PositionCount other)
        {
            int compareStart = Integer.compare(PosStartLower, other.PosStartLower);

            if(compareStart != 0)
                return compareStart;

            return Integer.compare(PosEndLower, other.PosEndLower);
        }

        public String toString()
        {
            return format("lower(%d-%d) upper(%d-%d) count(%d)", PosStartLower, PosStartUpper, PosEndLower, PosEndUpper, Count);
        }
    }
}
