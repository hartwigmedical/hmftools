package com.hartwig.hmftools.breakpointinspector.datamodel;

import java.util.Collection;

import com.hartwig.hmftools.breakpointinspector.Location;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.QueryInterval;

public class StructuralVariantResult {

    public Pair<Location, Location> Breakpoints;
    public Pair<Double, Double> AlleleFrequency = Pair.of(0.0, 0.0);
    public SampleStats TumorStats = new SampleStats();
    public SampleStats RefStats = new SampleStats();
    public Collection<String> Filters;
    public String FilterString = "";
    public QueryInterval[] QueryIntervals;
}
