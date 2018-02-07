package com.hartwig.hmftools.bachelor;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

class VariantModel {

    final VariantContext Context;
    final List<SnpEff> Annotations;
    final Set<String> dbSNP;

    List<SnpEff> SampleAnnotations;

    private static final Logger LOGGER = LogManager.getLogger(VariantModel.class);

    private VariantModel(final VariantContext ctx) {

        Context = ctx;
        dbSNP = Lists.newArrayList(ctx.getID().split(",")).stream().filter(s -> s.startsWith("rs")).collect(Collectors.toSet());
        Annotations = Arrays.stream(ctx.getAttributeAsString("ANN", "").split(","))
                .map(s -> Arrays.asList(s.split("\\|")))
                .map(SnpEff::parseAnnotation)
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
        
//        List<String> annotations = Lists.newArrayList(ctx.getAttributeAsString("ANN", "").split(","));
//
//        for(String annStr : annotations)
//        {
//            List<String> elements = Lists.newArrayList(annStr.split("\\|"));
//
//            if(elements.size() >= 6)
//            {
//                LOGGER.debug("size({}) e1={} e2={} e3={} e6={}", elements.size(), elements.get(0), elements.get(1), elements.get(2), elements.get(6));
//            }
//        }

        SampleAnnotations = Lists.newArrayList();
    }

    static VariantModel from(final VariantContext ctx) {
        return new VariantModel(ctx);
    }

    public void setSampleAnnotations(final List<String> alleleList) {
        SampleAnnotations.clear();

        SampleAnnotations = Annotations.stream()
                .filter(annotation -> alleleList.stream().anyMatch(allele -> allele.equals(annotation.getAllele())))
                .collect(Collectors.toList());
    }
}
