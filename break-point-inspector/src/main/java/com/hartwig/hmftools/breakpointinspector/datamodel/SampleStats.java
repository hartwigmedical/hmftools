package com.hartwig.hmftools.breakpointinspector.datamodel;

import static com.hartwig.hmftools.breakpointinspector.Util.prefixList;
import static com.hartwig.hmftools.breakpointinspector.Util.toStrings;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.breakpointinspector.clipping.Clipping;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class SampleStats {

    public final BreakpointStats bp1Stats = new BreakpointStats();
    public final BreakpointStats bp2Stats = new BreakpointStats();
    public final Clipping sampleClipping = new Clipping();
    public final List<Pair<SAMRecord, SAMRecord>> prEvidence = Lists.newArrayList();
    public final List<Pair<SAMRecord, SAMRecord>> srEvidence = Lists.newArrayList();

    @NotNull
    public static List<String> header() {
        final List<String> header = Lists.newArrayList();
        header.addAll(prefixList(BreakpointStats.HEADER, "BP1_"));
        header.addAll(prefixList(BreakpointStats.HEADER, "BP2_"));
        return header;
    }

    @NotNull
    public List<String> data() {
        final List<String> data = Lists.newArrayList();
        data.addAll(toStrings(data(bp1Stats)));
        data.addAll(toStrings(data(bp2Stats)));
        return data;
    }

    @NotNull
    private static List<Integer> data(@NotNull final BreakpointStats stats) {
        return Arrays.asList(stats.prOnlyNormal(), stats.prSrNormal(), stats.prOnlySupport(), stats.prSrSupport(), stats.srOnlySupport());
    }
}
