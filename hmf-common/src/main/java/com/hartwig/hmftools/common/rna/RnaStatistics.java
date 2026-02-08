package com.hartwig.hmftools.common.rna;

import java.util.List;

import org.immutables.value.Value;

@Value.Immutable
public abstract class RnaStatistics
{
    public abstract long totalFragments();
    public abstract long duplicateFragments();
    public abstract double splicedFragmentPerc();
    public abstract double unsplicedFragmentPerc();
    public abstract double altFragmentPerc();
    public abstract double chimericFragmentPerc();
    public abstract int splicedGeneCount();

    public abstract int readLength();

    // 5th, 50th and 95th percentile intronic fragment length
    public abstract double fragmentLength5thPercent();
    public abstract double fragmentLength50thPercent();
    public abstract double fragmentLength95thPercent();

    // proportion of fragments in 7 highly expressed genes
    public abstract double enrichedGenePercent();

    // Median GC (excluding 7 highly expressed genes)
    public abstract double medianGCRatio();

    public abstract double forwardStrandPercent();

    public abstract List<RnaQcFilter> qcStatus();
}
