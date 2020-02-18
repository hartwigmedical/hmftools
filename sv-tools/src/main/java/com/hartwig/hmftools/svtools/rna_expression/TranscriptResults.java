package com.hartwig.hmftools.svtools.rna_expression;

import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

import org.immutables.value.Value;

@Value.Immutable
public abstract class TranscriptResults
{
    public abstract TranscriptData trans();

    public abstract int exonsFound();
    public abstract int shortSupportingFragments();
    public abstract int shortUniqueFragments();
    public abstract int longSupportingFragments();
    public abstract int longUniqueFragments();
    public abstract int spliceJunctionsSupported();
    public abstract int spliceJunctionFragments();
    public abstract int spliceJunctionUniqueFragments();
    public abstract int uniqueSpliceJunctions();
    public abstract int uniqueSpliceJunctionsSupported();
    public abstract int exonicBases();
    public abstract int exonicBaseCoverage();
    public abstract int uniqueBases();
    public abstract int uniqueBaseCoverage();
    public abstract double uniqueBaseAvgDepth();
}
