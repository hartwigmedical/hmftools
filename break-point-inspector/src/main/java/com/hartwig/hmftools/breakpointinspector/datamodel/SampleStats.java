package com.hartwig.hmftools.breakpointinspector.datamodel;

import static com.hartwig.hmftools.breakpointinspector.Util.prefixList;
import static com.hartwig.hmftools.breakpointinspector.Util.toStrings;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.breakpointinspector.clipping.ClipInfo;
import com.hartwig.hmftools.breakpointinspector.clipping.Clipping;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class SampleStats {

    @NotNull
    private final BreakpointStats bp1Stats;
    @NotNull
    private final BreakpointStats bp2Stats;
    @NotNull
    private final List<Pair<SAMRecord, SAMRecord>> srEvidence;
    @NotNull
    private final List<Pair<SAMRecord, SAMRecord>> prEvidence;
    @NotNull
    private final Clipping sampleClipping = new Clipping();

    @NotNull
    static SampleStats emptyStats() {
        // KODU: This is fairly ugly, but needed to prevent NPE.
        return new SampleStats(new BreakpointStats(), new BreakpointStats(), Lists.newArrayList(), Lists.newArrayList());
    }
    
    public SampleStats(@NotNull final BreakpointStats bp1Stats, @NotNull final BreakpointStats bp2Stats,
            @NotNull final List<Pair<SAMRecord, SAMRecord>> srEvidence, @NotNull final List<Pair<SAMRecord, SAMRecord>> prEvidence) {
        this.bp1Stats = bp1Stats;
        this.bp2Stats = bp2Stats;
        this.srEvidence = srEvidence;
        this.prEvidence = prEvidence;
    }

    public void addSampleClipping(@NotNull final ClipInfo clipInfo) {
        sampleClipping.add(clipInfo);
    }

    @NotNull
    public BreakpointStats bp1Stats() {
        return bp1Stats;
    }

    @NotNull
    public BreakpointStats bp2Stats() {
        return bp2Stats;
    }

    @NotNull
    public List<Pair<SAMRecord, SAMRecord>> srEvidence() {
        return srEvidence;
    }

    @NotNull
    public List<Pair<SAMRecord, SAMRecord>> prEvidence() {
        return prEvidence;
    }

    @NotNull
    public Clipping sampleClipping() {
        return sampleClipping;
    }

    @NotNull
    public static List<String> header() {
        final List<String> header = Lists.newArrayList();
        header.addAll(prefixList(BreakpointStats.HEADER, "BP1_"));
        header.addAll(prefixList(BreakpointStats.HEADER, "BP2_"));
        return header;
    }

    @NotNull
    public List<String> statsData() {
        final List<String> data = Lists.newArrayList();
        data.addAll(toStrings(statsData(bp1Stats)));
        data.addAll(toStrings(statsData(bp2Stats)));
        return data;
    }

    @NotNull
    private static List<Integer> statsData(@NotNull final BreakpointStats stats) {
        return Arrays.asList(stats.prOnlyNormal(), stats.prSrNormal(), stats.prOnlySupport(), stats.prSrSupport(), stats.srOnlySupport());
    }
}
