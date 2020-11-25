package com.hartwig.hmftools.common.germline;

import com.hartwig.hmftools.common.variant.CodingEffect;

import org.immutables.value.Value;

@Value.Immutable
public abstract class GermlineVariant
{
    // database fields
    public abstract String chromosome();
    public abstract long position();
    public abstract FilterType filter();
    public abstract boolean reported();
    public abstract PathogenicType pathogenic();
    public abstract String type();
    public abstract String ref();
    public abstract String alts();
    public abstract String gene();
    public abstract String transcriptId();
    public abstract String effects();
    public abstract CodingEffect codingEffect();
    public abstract String microhomology();
    public abstract String repeatSequence();
    public abstract int repeatCount();
    public abstract int alleleReadCount();
    public abstract int totalReadCount();
    public abstract double adjustedVaf();
    public abstract double adjustedCopyNumber();
    public abstract String trinucleotideContext();
    public abstract String hgvsProtein();
    public abstract String hgvsCoding();
    public abstract boolean biallelic();
    public abstract double minorAlleleJcn();
    public abstract String refStatus();
    public abstract String clinvarInfo();

}
