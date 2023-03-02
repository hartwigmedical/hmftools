package com.hartwig.hmftools.datamodel.hla;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class LilacAllele {

    public abstract String allele();
    public abstract int refFragments();
    public abstract int refUnique();
    public abstract int refShared();
    public abstract int refWild();
    public abstract int tumorFragments();
    public abstract int tumorUnique();
    public abstract int tumorShared();
    public abstract int tumorWild();
    public abstract int rnaFragments();
    public abstract int rnaUnique();
    public abstract int rnaShared();
    public abstract int rnaWild();

    public abstract double tumorCopyNumber();

    public abstract double somaticMissense();
    public abstract double somaticNonsenseOrFrameshift();
    public abstract double somaticSplice();
    public abstract double somaticSynonymous();
    public abstract double somaticInframeIndel();

    public double somaticVariantCount() {
        return somaticMissense() + somaticNonsenseOrFrameshift() + somaticSplice() + somaticSynonymous() + somaticInframeIndel();
    }
}
