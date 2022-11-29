package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.Variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableVariant implements Variant {

    @NotNull
    public abstract ReportableVariantSource source();

    @NotNull
    @Override
    public abstract String gene();

    @NotNull
    public abstract String transcript();

    public abstract boolean isCanonical();

    @NotNull
    @Override
    public abstract String chromosome();

    @Override
    public abstract int position();

    @NotNull
    @Value.Derived
    public String gDNA()
    {
        return chromosome() + ":" + position();
    }

    @NotNull
    @Override
    public abstract String ref();

    @NotNull
    @Override
    public abstract String alt();

    @NotNull
    @Override
    public abstract String canonicalTranscript();

    @NotNull
    @Override
    public abstract String canonicalEffect();

    @NotNull
    @Override
    public abstract CodingEffect canonicalCodingEffect();

    @NotNull
    @Override
    public abstract String canonicalHgvsCodingImpact();

    @NotNull
    @Override
    public abstract String canonicalHgvsProteinImpact();

    @NotNull
    public abstract String otherReportedEffects();

    @Override
    public abstract int totalReadCount();

    @Override
    public abstract int alleleReadCount();

    @Nullable
    public abstract Integer rnaAlleleReadCount();

    @Nullable
    public abstract Integer rnaTotalReadCount();

    public abstract double totalCopyNumber();

    public abstract double alleleCopyNumber();

    public abstract double minorAlleleCopyNumber();

    @NotNull
    @Value.Derived
    public String tVAF()
    {
        double vaf = alleleCopyNumber() / totalCopyNumber();
        return DataUtil.formatPercentage(100 * Math.max(0, Math.min(1, vaf)));
    }

    @NotNull
    public abstract Hotspot hotspot();

    public abstract double clonalLikelihood();

    public abstract double driverLikelihood();

    @NotNull
    @Value.Derived
    public DriverInterpretation driverLikelihoodInterpretation()
    {
        return DriverInterpretation.interpret(driverLikelihood());
    }

    @Nullable
    public abstract Boolean biallelic();

    @NotNull
    public abstract GenotypeStatus genotypeStatus();

    @Nullable
    public abstract Integer localPhaseSet();
}