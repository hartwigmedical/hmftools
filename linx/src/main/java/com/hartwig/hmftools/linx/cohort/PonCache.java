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
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.cohort.PonBuilder.FLD_COUNT;
import static com.hartwig.hmftools.linx.cohort.PonBuilder.formChrOrientKey;
import static com.hartwig.hmftools.linx.cohort.PonMatchType.EXACT;
import static com.hartwig.hmftools.linx.cohort.PonMatchType.MARGIN;
import static com.hartwig.hmftools.linx.cohort.PonMatchType.NONE;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;

public class PonCache
{
    private final Map<String,List<SvPonEntry>> mSvPositionsMap;
    private final Map<String,List<SglPonEntry>> mSglPositionsMap;
    private boolean mEnabledForAnnotation;

    public PonCache()
    {
        mSvPositionsMap = Maps.newHashMap();
        mSglPositionsMap = Maps.newHashMap();
        mEnabledForAnnotation = false;
    }

    public boolean hasEntries() { return mEnabledForAnnotation; }

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

                int count = Integer.parseInt(values[countIndex]);

                if(values[chrEndIndex].isEmpty())
                {
                    String chrOrientKey = formChrOrientKey(chrStart, orientStart, null, null);
                    List<SglPonEntry> positions = mSglPositionsMap.get(chrOrientKey);

                    if(positions == null)
                    {
                        positions = Lists.newArrayList();
                        mSglPositionsMap.put(chrOrientKey, positions);
                    }

                    positions.add(new SglPonEntry(posStart, count));
                }
                else
                {
                    String chrEnd = values[chrEndIndex];
                    int posEnd = Integer.parseInt(values[posEndIndex]);
                    Orientation orientEnd = Orientation.fromByteStr(values[orientEndIndex]);

                    String chrOrientKey = formChrOrientKey(chrStart, orientStart, chrEnd, orientEnd);
                    List<SvPonEntry> svPonEntries = mSvPositionsMap.get(chrOrientKey);

                    if(svPonEntries == null)
                    {
                        svPonEntries = Lists.newArrayList();
                        mSvPositionsMap.put(chrOrientKey, svPonEntries);
                    }

                    svPonEntries.add(new SvPonEntry(posStart, posEnd, count));
                }
            }

            LNX_LOGGER.debug("loaded {} PON entries", entries);

            // sort each map by position
            for(List<SglPonEntry> ponEntries : mSglPositionsMap.values())
            {
                Collections.sort(ponEntries);
            }

            for(List<SvPonEntry> svPonEntries : mSvPositionsMap.values())
            {
                Collections.sort(svPonEntries);
            }

            mEnabledForAnnotation = true;
        }
        catch(IOException exception)
        {
            LNX_LOGGER.error("failed to read PON file({})", filename, exception.toString());
        }
    }

    public PonMatch matchesSglEntry(final String chrStart, final int posStart, final Orientation orientStart, int margin)
    {
        String chrOrientKey = formChrOrientKey(chrStart, orientStart, null, null);
        List<SglPonEntry> ponEntries = mSglPositionsMap.get(chrOrientKey);

        if(ponEntries == null)
            return PonMatch.NONE;

        PonMatch ponMatch = PonMatch.NONE;

        for(SglPonEntry ponEntry : ponEntries)
        {
            if(posStart > ponEntry.Position + margin)
                continue;

            if(posStart < ponEntry.Position - margin)
                break;

            if(ponEntry.Position == posStart)
                return new PonMatch(EXACT, ponEntry.Count);

            if(abs(ponEntry.Position - posStart) <= margin)
                ponMatch = new PonMatch(MARGIN, ponEntry.Count);
        }

        return ponMatch;
    }

    public PonMatch matchesSvEntry(
            final String chrStart, final int posStart, final Orientation orientStart,
            final String chrEnd, final int posEnd, final Orientation orientEnd, int margin)
    {
        String chrOrientKey = formChrOrientKey(chrStart, orientStart, chrEnd, orientEnd);
        List<SvPonEntry> svPonEntries = mSvPositionsMap.get(chrOrientKey);

        if(svPonEntries == null)
            return PonMatch.NONE;

        PonMatch ponMatch = PonMatch.NONE;

        for(SvPonEntry ponEntry : svPonEntries)
        {
            if(posStart > ponEntry.PosStart + margin)
                continue;

            if(posStart < ponEntry.PosStart - margin)
                break;

            if(ponEntry.PosStart == posStart && ponEntry.PosEnd == posEnd)
                return new PonMatch(EXACT, ponEntry.Count);

            if(abs(ponEntry.PosStart - posStart) <= margin && abs(ponEntry.PosEnd - posEnd) <= margin)
                ponMatch = new PonMatch(MARGIN, ponEntry.Count);
        }

        return ponMatch;
    }

    private class SglPonEntry implements Comparable<SglPonEntry>
    {
        public final int Position;
        public final int Count;

        public SglPonEntry(final int position, final int count)
        {
            Position = position;
            Count = count;
        }

        @Override
        public int compareTo(final SglPonEntry other)
        {
            return Integer.compare(Position, other.Position);
        }

        public String toString() { return format("%d_%d", Position, Count); }
    }

    private class SvPonEntry implements Comparable<SvPonEntry>
    {
        public final int PosStart;
        public final int PosEnd;
        public final int Count;

        public SvPonEntry(final int posStart, final int posEnd, final int count)
        {
            PosStart = posStart;
            PosEnd = posEnd;
            Count = count;
        }

        @Override
        public int compareTo(final SvPonEntry other)
        {
            int compareStart = Integer.compare(PosStart, other.PosStart);

            if(compareStart != 0)
                return compareStart;

            return Integer.compare(PosEnd, other.PosEnd);
        }

        public String toString() { return format("%d_%d", PosStart, PosEnd); }
    }
}
