package com.hartwig.hmftools.common.variant;

import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlProperty;
import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.utils.DataUtil;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableVariant implements Variant {

    @NotNull
    public abstract ReportableVariantSource source();

    @JacksonXmlProperty(localName = "name")
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

    @JacksonXmlProperty(localName = "pos")
    @NotNull
    @Value.Derived
    public String gDNA() {
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

    @JacksonXmlProperty(localName = "var")
    @NotNull
    @Override
    public abstract String canonicalHgvsCodingImpact();

    @JacksonXmlProperty(localName = "prot")
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

    @JacksonXmlProperty(localName = "tvaf")
    @NotNull
    @Value.Derived
    public String tVAF() {
        double vaf = alleleCopyNumber() / totalCopyNumber();
        return DataUtil.formatPercentage(100 * Math.max(0, Math.min(1, vaf)));
    }

    @JacksonXmlProperty(localName = "hotsp")
    @NotNull
    public abstract Hotspot hotspot();

    public abstract double clonalLikelihood();

    public abstract double driverLikelihood();

    @JacksonXmlProperty(localName = "driver")
    @NotNull
    @Value.Derived
    public DriverInterpretation driverLikelihoodInterpretation() {
        return DriverInterpretation.interpret(driverLikelihood());
    }

    @JacksonXmlProperty(localName = "biallic")
    @Nullable
    public abstract Boolean biallelic();

    @NotNull
    public abstract GenotypeStatus genotypeStatus();

    @Nullable
    public abstract Integer localPhaseSet();
}
