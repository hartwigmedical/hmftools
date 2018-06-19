package com.hartwig.hmftools.breakpointinspector.datamodel;

import java.util.Collection;

import com.hartwig.hmftools.breakpointinspector.Location;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.QueryInterval;

public class StructuralVariantResult {

    public Pair<Location, Location> breakpoints;
    public Pair<Double, Double> alleleFrequency = Pair.of(0.0, 0.0);
    public SampleStats tumorStats = new SampleStats();
    public SampleStats refStats = new SampleStats();
    public Collection<String> filters;
    public String filterString = "";
    public QueryInterval[] queryIntervals;
}
