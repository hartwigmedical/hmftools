package com.hartwig.hmftools.linx.cohort;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHR_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHR_START;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENT_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENT_START;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.cohort.PonMatchType.NONE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;

public class PonCache
{
    // used for building the PON
    // keyed by chrStart_orientStart_chrEnd_orientEnd,then posStart_posEnd for SVs
    private final Map<String,Map<String,Integer>> mSvPositionCountsMap;

    // keyed by chrStart_orientStart, then posStart for SGLs
    private final Map<String,Map<Integer,Integer>> mSglPositionCountsMap;

    // used for applying the PON
    private final Map<String, List<PositionPair>> mSvPositionsMap;
    private final Map<String,List<Integer>> mSglPositionsMap;
    private boolean mEnabledForAnnotation;

    public PonCache()
    {
        mSvPositionCountsMap = Maps.newHashMap();
        mSglPositionCountsMap = Maps.newHashMap();

        mSvPositionsMap = Maps.newHashMap();
        mSglPositionsMap = Maps.newHashMap();
        mEnabledForAnnotation = false;
    }

    // public Map<String,Map<String,Integer>> chrLocationCountsMap() { return mSvPositionCountsMap; }
    public boolean hasEntries() { return mEnabledForAnnotation; }

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

    private static String PAIR_DELIM = "_";

    private static String formChrOrientKey(
            final String chrStart, final Orientation orientStart, @Nullable final String chrEnd, @Nullable final Orientation orientEnd)
    {
        if(chrEnd == null)
            return format("%s_%d", chrStart, orientStart.asByte());
        else
            return format("%s_%d_%s_%d", chrStart, orientStart.asByte(), chrEnd, orientEnd.asByte());
    }

    private static String formPositionPair(final int posStart, final Integer posEnd)
    {
        return format("%d_%d", posStart, posEnd);
    }

    private static final String FLD_COUNT = "Count";

    public void writePonCache(final CohortConfig config)
    {
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

                for(Map.Entry<String,Integer> locationEntry : locationMap.entrySet())
                {
                    int count = locationEntry.getValue();

                    if(count == 1)
                        continue;

                    String positionKey = locationEntry.getKey();

                    String[] positionParts = positionKey.split(PAIR_DELIM, 2);

                    sj = new StringJoiner(TSV_DELIM);

                    sj.add(chrStart).add(positionParts[0]).add(orientStart);
                    sj.add(chrEnd).add(positionParts[1]).add(orientEnd);

                    sj.add(String.valueOf(count));

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

                for(Map.Entry<Integer,Integer> positionEntry : locationMap.entrySet())
                {
                    int posStart = positionEntry.getKey();
                    int count = positionEntry.getValue();

                    if(count == 1)
                        continue;

                    sj = new StringJoiner(TSV_DELIM);

                    sj.add(chrStart).add(String.valueOf(posStart)).add(orientStart);
                    sj.add("").add("").add("");

                    sj.add(String.valueOf(count));

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

    public void loadPonFile(final String filename)
    {
        if(filename == null || filename.isEmpty())
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String header = fileReader.readLine();
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

            int chrStartIndex = fieldsIndexMap.get(FLD_CHR_START);
            int chrEndIndex = fieldsIndexMap.get(FLD_CHR_END);
            int posStartIndex = fieldsIndexMap.get(FLD_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_POS_END);
            int orientStartIndex = fieldsIndexMap.get(FLD_ORIENT_START);
            int orientEndIndex = fieldsIndexMap.get(FLD_ORIENT_END);
            int countIndex = fieldsIndexMap.get(FLD_COUNT);

            String line;

            int entries = 0;

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(TSV_DELIM, -1);

                ++entries;

                String chrStart = values[chrStartIndex];
                int posStart = Integer.parseInt(values[posStartIndex]);
                Orientation orientStart = Orientation.fromByteStr(values[orientStartIndex]);

                // int count = Integer.parseInt(values[countIndex]);

                if(values[chrEndIndex].isEmpty())
                {
                    String chrOrientKey = formChrOrientKey(chrStart, orientStart, null, null);
                    List<Integer> positions = mSglPositionsMap.get(chrOrientKey);

                    if(positions == null)
                    {
                        positions = Lists.newArrayList();
                        mSglPositionsMap.put(chrOrientKey, positions);
                    }

                    positions.add(posStart);
                }
                else
                {
                    String chrEnd = values[chrEndIndex];
                    int posEnd = Integer.parseInt(values[posEndIndex]);
                    Orientation orientEnd = Orientation.fromByteStr(values[orientEndIndex]);

                    String chrOrientKey = formChrOrientKey(chrStart, orientStart, chrEnd, orientEnd);
                    List<PositionPair> positionPairs = mSvPositionsMap.get(chrOrientKey);

                    if(positionPairs == null)
                    {
                        positionPairs = Lists.newArrayList();
                        mSvPositionsMap.put(chrOrientKey, positionPairs);
                    }

                    positionPairs.add(new PositionPair(posStart, posEnd));
                }
            }

            LNX_LOGGER.debug("loaded {} PON entries", entries);

            // sort each map by position
            for(List<Integer> positions : mSglPositionsMap.values())
            {
                Collections.sort(positions);
            }

            for(List<PositionPair> positionPairs : mSvPositionsMap.values())
            {
                Collections.sort(positionPairs);
            }

            mEnabledForAnnotation = true;
        }
        catch(IOException exception)
        {
            LNX_LOGGER.error("failed to read PON file({})", filename, exception.toString());
        }
    }

    public PonMatchType matchesSglEntry(final String chrStart, final int posStart, final Orientation orientStart, int margin)
    {
        String chrOrientKey = formChrOrientKey(chrStart, orientStart, null, null);
        List<Integer> positions = mSglPositionsMap.get(chrOrientKey);

        if(positions == null)
            return NONE;

        PonMatchType ponMatchType = NONE;

        for(int position : positions)
        {
            if(posStart > position + margin)
                continue;

            if(posStart < position - margin)
                break;

            if(position == posStart)
                return PonMatchType.EXACT;

            if(abs(position - posStart) <= margin)
                ponMatchType = PonMatchType.MARGIN;
        }

        return ponMatchType;
    }

    public PonMatchType matchesSvEntry(
            final String chrStart, final int posStart, final Orientation orientStart,
            final String chrEnd, final int posEnd, final Orientation orientEnd, int margin)
    {
        String chrOrientKey = formChrOrientKey(chrStart, orientStart, chrEnd, orientEnd);
        List<PositionPair> positionPairs = mSvPositionsMap.get(chrOrientKey);

        if(positionPairs == null)
            return NONE;

        PonMatchType ponMatchType = NONE;

        for(PositionPair positionPair : positionPairs)
        {
            if(posStart > positionPair.PosStart + margin)
                continue;

            if(posStart < positionPair.PosStart - margin)
                break;

            if(positionPair.PosStart == posStart && positionPair.PosEnd == posEnd)
                return PonMatchType.EXACT;

            if(abs(positionPair.PosStart - posStart) <= margin && abs(positionPair.PosEnd - posEnd) <= margin)
                ponMatchType = PonMatchType.MARGIN;
        }

        return ponMatchType;
    }

    private class PositionPair implements Comparable<PositionPair>
    {
        public final int PosStart;
        public final int PosEnd;

        public PositionPair(final int posStart, final int posEnd)
        {
            PosStart = posStart;
            PosEnd = posEnd;
        }

        @Override
        public int compareTo(final PositionPair other)
        {
            int compareStart = Integer.compare(PosStart, other.PosStart);

            if(compareStart != 0)
                return compareStart;

            return Integer.compare(PosEnd, other.PosEnd);
        }

        public String toString() { return format("%d_%d", PosStart, PosEnd); }
    }
}
