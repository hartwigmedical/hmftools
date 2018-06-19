package com.hartwig.hmftools.breakpointinspector.datamodel;

import static com.hartwig.hmftools.breakpointinspector.Util.prefixList;
import static com.hartwig.hmftools.breakpointinspector.Util.toStrings;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.breakpointinspector.clipping.Clipping;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class SampleStats {

    public final BreakpointStats BP1_Stats = new BreakpointStats();
    public final BreakpointStats BP2_Stats = new BreakpointStats();
    public final Clipping Sample_Clipping = new Clipping();
    public final List<Pair<SAMRecord, SAMRecord>> PR_Evidence = Lists.newArrayList();
    public final List<Pair<SAMRecord, SAMRecord>> SR_Evidence = Lists.newArrayList();

    @NotNull
    public static List<String> GetHeader() {
        final List<String> header = Lists.newArrayList();
        header.addAll(prefixList(BreakpointStats.GetHeader(), "BP1_"));
        header.addAll(prefixList(BreakpointStats.GetHeader(), "BP2_"));
        return header;
    }

    @NotNull
    public List<String> GetData() {
        final List<String> data = Lists.newArrayList();
        data.addAll(toStrings(BP1_Stats.GetData()));
        data.addAll(toStrings(BP2_Stats.GetData()));
        return data;
    }
}
