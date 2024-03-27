package com.hartwig.hmftools.cobalt.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;

public class FragmentGcMap
{
    private final Map<String,FragmentGcCounts> mMap;

    public FragmentGcMap()
    {
        mMap = Maps.newHashMap();
    }

    public void add(final int fragmentLength, final double gcPercent, final int duplicateCount)
    {
        String key = FragmentGcCounts.formKey(fragmentLength, duplicateCount, gcPercent);

        FragmentGcCounts counts = mMap.get(key);

        if(counts == null)
        {
            counts = new FragmentGcCounts(fragmentLength, duplicateCount, gcPercent);
            mMap.put(key, counts);
        }

        ++counts.Count;
    }

    public void merge(final FragmentGcMap other)
    {
        for(Map.Entry<String,FragmentGcCounts> entry : other.mMap.entrySet())
        {
            String key = entry.getKey();
            FragmentGcCounts otherCounts = entry.getValue();

            FragmentGcCounts fragGcCounts = mMap.get(key);

            if(fragGcCounts == null)
            {
                mMap.put(key, otherCounts);
            }
            else
            {
                fragGcCounts.Count += otherCounts.Count;
            }
        }
    }

    public void writeMap(final BufferedWriter writer, final String regionStr) throws IOException
    {
        if(mMap.isEmpty())
            return;

        List<FragmentGcCounts> fragmentGcCounts = mMap.values().stream().collect(Collectors.toList());
        Collections.sort(fragmentGcCounts);

        for(FragmentGcCounts fragGcCount : fragmentGcCounts)
        {
            StringJoiner sb = new StringJoiner(TSV_DELIM);
            sb.add(regionStr);
            sb.add(String.valueOf(fragGcCount.FragmentLength));
            sb.add(String.valueOf(fragGcCount.DuplicateCount));
            sb.add(format("%.2f", fragGcCount.GcPercent));
            sb.add(String.valueOf(fragGcCount.Count));
            writer.write(sb.toString());
            writer.newLine();
        }
    }
}
