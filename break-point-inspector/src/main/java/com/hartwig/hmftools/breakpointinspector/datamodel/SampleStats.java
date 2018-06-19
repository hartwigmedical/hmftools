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

    public final BreakpointStats bp1Stats = new BreakpointStats();
    public final BreakpointStats bp2Stats = new BreakpointStats();
    public final Clipping sampleClipping = new Clipping();
    public final List<Pair<SAMRecord, SAMRecord>> prEvidence = Lists.newArrayList();
    public final List<Pair<SAMRecord, SAMRecord>> srEvidence = Lists.newArrayList();

    @NotNull
    public static List<String> header() {
        final List<String> header = Lists.newArrayList();
        header.addAll(prefixList(BreakpointStats.header(), "BP1_"));
        header.addAll(prefixList(BreakpointStats.header(), "BP2_"));
        return header;
    }

    @NotNull
    public List<String> data() {
        final List<String> data = Lists.newArrayList();
        data.addAll(toStrings(bp1Stats.data()));
        data.addAll(toStrings(bp2Stats.data()));
        return data;
    }
}
