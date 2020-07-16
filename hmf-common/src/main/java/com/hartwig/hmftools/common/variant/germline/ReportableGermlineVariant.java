package com.hartwig.hmftools.common.variant.germline;

import com.hartwig.hmftools.common.variant.CodingEffect;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableGermlineVariant
{
    public abstract String passFilter();
    public abstract String gene();
    public abstract String ref();
    public abstract String alt();
    public abstract CodingEffect codingEffect();
    public abstract String chromosome();
    public abstract long position();
    public abstract String hgvsProtein();
    public abstract String hgvsCoding();
    public abstract int totalReadCount();
    public abstract int alleleReadCount();
    public abstract double adjustedVaf();
    public abstract double adjustedCopyNumber();
    public abstract boolean biallelic();

}
