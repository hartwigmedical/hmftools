package com.hartwig.hmftools.bachelor;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.snpeff.VariantAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.VariantAnnotationFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantModel {

    private final VariantContext context;
    private final List<VariantAnnotation> annotations;
    private final Set<String> dbSNP;
    private final List<VariantAnnotation> sampleAnnotations;

    private static final Logger LOGGER = LogManager.getLogger(VariantModel.class);

    public VariantModel(final String sample, final VariantContext ctx) {

        context = ctx;
        dbSNP = Lists.newArrayList(ctx.getID().split(",")).stream().filter(s -> s.startsWith("rs")).collect(Collectors.toSet());
        annotations = VariantAnnotationFactory.fromContext(ctx);

        final List<String> alleleList =
                ctx.getGenotype(sample).getAlleles().stream().map(Allele::getBaseString).collect(Collectors.toList());
        sampleAnnotations = annotations.stream()
                .filter(annotation -> alleleList.stream().anyMatch(allele -> allele.equals(annotation.allele())))
                .collect(Collectors.toList());

        for (String allele : alleleList) {
            LOGGER.debug("checking allele({}):", allele);
        }
        LOGGER.debug("annotation alleleCount(reduced={} orig={}) v listCount({}):",
                sampleAnnotations().size(),
                annotations().size(),
                alleleList.size());
    }

    public VariantContext context() {
        return context;
    }

    public List<VariantAnnotation> annotations() {
        return annotations;
    }

    public Set<String> dbSNP() {
        return dbSNP;
    }

    public List<VariantAnnotation> sampleAnnotations() {
        return sampleAnnotations;
    }
}
