package com.hartwig.hmftools.common.rna;

import org.immutables.value.Value;

@Value.Immutable
public abstract class TranscriptExpression
{
    public abstract String transcriptName();
    public abstract String geneName();
    public abstract double tpm();
}
