package com.hartwig.hmftools.redux.merge;

import static java.lang.String.format;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import htsjdk.samtools.QueryInterval;

public class SequenceInfo
{
    public final int Id;
    public final String BamFile;
    public final List<QueryInterval> Intervals;

    public SequenceInfo(final int id, final String bamFile)
    {
        Id = id;
        BamFile = bamFile;
        Intervals = Lists.newArrayList();
    }

    public QueryInterval[] asArray()
    {
        QueryInterval[] intervals = new QueryInterval[Intervals.size()];

        for(int i = 0; i < Intervals.size(); ++i)
        {
            intervals[i] = Intervals.get(i);
        }

        return intervals;
    }

    public long intervalLength() { return Intervals.stream().mapToLong(x -> x.end-  x.start + 1).sum(); }

    public String toString()
    {
        StringJoiner sj = new StringJoiner(", ");

        for(QueryInterval interval : Intervals)
        {
            sj.add(format("%d:%d-%d", interval.referenceIndex, interval.start, interval.end));
        }

        return format("%d: length(%d) intervals: %s", Id, intervalLength(), sj);
    }
}
