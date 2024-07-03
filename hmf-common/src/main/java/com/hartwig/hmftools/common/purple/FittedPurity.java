package com.hartwig.hmftools.common.purple;

import java.util.Comparator;
import java.util.List;
import java.util.TreeSet;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;

import org.immutables.value.Value;

@Value.Immutable
public abstract class FittedPurity implements Comparable<FittedPurity>
{
    public abstract double purity();

    public abstract double normFactor();

    public abstract double ploidy();

    public abstract double score();

    public abstract double diploidProportion();

    public abstract double somaticPenalty();

    @Override
    public int compareTo(final FittedPurity o) {
        return Double.compare(score(), o.score());
    }
}
