package com.hartwig.hmftools.bachelor;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.snpeff.VariantAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.VariantAnnotationFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantModel {

    private static final Logger LOGGER = LogManager.getLogger(VariantModel.class);

    private final VariantContext context;
    private final Set<String> dbSNP;
    private final List<VariantAnnotation> sampleAnnotations;

    public VariantModel(@NotNull String sample, @NotNull VariantContext context) {
        this.context = context;
        dbSNP = Lists.newArrayList(context.getID().split(",")).stream().filter(s -> s.startsWith("rs")).collect(Collectors.toSet());
        final List<VariantAnnotation> annotations = VariantAnnotationFactory.fromContext(context);

        final List<String> alleleList =
                context.getGenotype(sample).getAlleles().stream().map(Allele::getBaseString).collect(Collectors.toList());
        sampleAnnotations = annotations.stream()
                .filter(annotation -> alleleList.stream().anyMatch(allele -> allele.equals(annotation.allele())))
                .collect(Collectors.toList());

//        LOGGER.debug("annotation alleleCount(reduced={} orig={}) v listCount({}):", sampleAnnotations.size(), annotations.size(),
//                alleleList.size());
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
    public List<VariantAnnotation> sampleAnnotations() {
        return sampleAnnotations;
    }
}
