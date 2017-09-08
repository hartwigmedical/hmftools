package com.hartwig.hmftools.breakpointinspector;

import static com.hartwig.hmftools.breakpointinspector.Util.prefixList;
import static com.hartwig.hmftools.breakpointinspector.Util.toStrings;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.breakpointinspector.clipping.Clipping;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;

class BreakpointStats {
    int PR_Only_Normal = 0;
    int PR_SR_Normal = 0;
    int PR_Only_Support = 0;
    int PR_SR_Support = 0;
    int SR_Only_Support = 0;

    static List<String> GetHeader() {
        return Arrays.asList("PR_ONLY_NORMAL", "PR_SR_NORMAL", "PR_ONLY_SUPPORT", "PR_SR_SUPPORT", "SR_ONLY_SUPPORT");
    }

    List<Integer> GetData() {
        return Arrays.asList(PR_Only_Normal, PR_SR_Normal, PR_Only_Support, PR_SR_Support, SR_Only_Support);
    }
}

class SampleStats {
    BreakpointStats BP1_Stats = new BreakpointStats();
    BreakpointStats BP2_Stats = new BreakpointStats();
    Clipping Sample_Clipping = new Clipping();
    List<Pair<SAMRecord, SAMRecord>> PR_Evidence = Lists.newArrayList();

    static List<String> GetHeader() {
        final List<String> header = Lists.newArrayList();
        header.addAll(prefixList(BreakpointStats.GetHeader(), "BP1_"));
        header.addAll(prefixList(BreakpointStats.GetHeader(), "BP2_"));
        return header;
    }

    List<String> GetData() {
        final List<String> data = Lists.newArrayList();
        data.addAll(toStrings(BP1_Stats.GetData()));
        data.addAll(toStrings(BP2_Stats.GetData()));
        return data;
    }
}

class StructuralVariantResult {
    Pair<Location, Location> Breakpoints;
    Pair<Double, Double> AlleleFrequency = Pair.of(0.0, 0.0);
    SampleStats TumorStats = new SampleStats();
    SampleStats RefStats = new SampleStats();
    Collection<String> Filters;
    String FilterString = "";
    QueryInterval[] QueryIntervals;
}