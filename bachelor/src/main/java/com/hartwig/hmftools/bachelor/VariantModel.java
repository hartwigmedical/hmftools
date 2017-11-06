package com.hartwig.hmftools.bachelor;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import htsjdk.variant.variantcontext.VariantContext;

class VariantModel {

    final VariantContext Context;
    final List<SnpEff> Annotations;
    final Set<String> dbSNP;

    private VariantModel(final VariantContext ctx) {
        Context = ctx;
        dbSNP = Lists.newArrayList(ctx.getID().split(",")).stream().filter(s -> s.startsWith("rs")).collect(Collectors.toSet());
        Annotations = Arrays.stream(ctx.getAttributeAsString("ANN", "").split(","))
                .map(s -> Arrays.asList(s.split("\\|")))
                .map(SnpEff::parseAnnotation)
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
    }

    static VariantModel from(final VariantContext ctx) {
        return new VariantModel(ctx);
    }

}
