package com.hartwig.hmftools.bachelor.types;

import com.hartwig.hmftools.common.variant.CodingEffect;

import org.immutables.value.Value;

@Value.Immutable
public abstract class GermlineVariant
{
    // database fields
    public abstract String chromosome();
    public abstract long position();
    public abstract FilterType filter();
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

    // other fields
    public abstract String program();
    public abstract String variantId();
    public abstract String annotations();
    public abstract int phredScore();
    public abstract boolean isHomozygous();
    public abstract String matchType();
    public abstract String codonInfo();
    public abstract boolean clinvarMatch();
    public abstract String clinvarSignificance();
    public abstract String clinvarSignificanceInfo();

}
