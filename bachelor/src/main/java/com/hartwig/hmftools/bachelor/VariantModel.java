package com.hartwig.hmftools.bachelor;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantModel {

    private final VariantContext context;
    private final Set<String> dbSNP;
    private final List<SnpEffAnnotation> sampleAnnotations;
    private final List<String> rawAnnotations;

    public VariantModel(@NotNull String sample, @NotNull VariantContext context) {
        this.context = context;
        dbSNP = Lists.newArrayList(context.getID().split(",")).stream().filter(s -> s.startsWith("rs")).collect(Collectors.toSet());
        final List<SnpEffAnnotation> annotations = SnpEffAnnotationFactory.fromContext(context);

        rawAnnotations = SnpEffAnnotationFactory.rawAnnotations(context);

        final List<String> alleleList =
                context.getGenotype(sample).getAlleles().stream().map(Allele::getBaseString).collect(Collectors.toList());
        sampleAnnotations = annotations.stream()
                .filter(annotation -> alleleList.stream().anyMatch(allele -> allele.equals(annotation.allele())))
                .collect(Collectors.toList());
    }

    @NotNull
    public VariantContext context() {
        return context;
    }

    @NotNull
    public Set<String> dbSNP() {
        return dbSNP;
    }

    @NotNull
    public List<SnpEffAnnotation> sampleAnnotations() {
        return sampleAnnotations;
    }

    @NotNull
    public List<String> rawAnnotations() {
        return rawAnnotations;
    }
}
