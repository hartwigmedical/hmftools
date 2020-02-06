package com.hartwig.hmftools.svtools.rna_expression;

import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

import org.immutables.value.Value;

@Value.Immutable
public abstract class TranscriptResults
{
    public abstract TranscriptData trans();

    public abstract int exonsFound();
    public abstract int supportingReads();
    public abstract int spliceJunctionsSupported();
    public abstract int spliceJunctionReads();
    public abstract double avgReadDepth();
    public abstract long codingLengthTotal();
    public abstract long totalReadCoverage();
    public abstract int intronicReads();
    public abstract int uniqueIntronicReads();

}
