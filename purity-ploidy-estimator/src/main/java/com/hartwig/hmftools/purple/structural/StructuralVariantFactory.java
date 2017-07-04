package com.hartwig.hmftools.purple.structural;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

class StructuralVariantFactory {

    private static String TYPE = "SVTYPE";
    private static String MATE_ID = "MATEID";
    private static String INS_SEQ = "SVINSSEQ";
    private static String HOM_SEQ = "HOMSEQ";

    private final Map<String, VariantContext> unmatched = new HashMap<>();
    private final List<StructuralVariant> results = Lists.newArrayList();
    private final PassingVariantFilter filter = new PassingVariantFilter();
    private final Set<String> samples = Sets.newHashSet();

    void addVariantContext(VariantContext context) {
        if (filter.test(context)) {

            samples.addAll(context.getSampleNames());
            final StructuralVariantType type = type(context);
            if (type.equals(StructuralVariantType.BND)) {
                final String mate = (String) context.getAttribute(MATE_ID);
                if (unmatched.containsKey(mate)) {
                    results.add(create(unmatched.remove(mate), context));
                } else {
                    unmatched.put(context.getID(), context);
                }
            } else {
                results.add(create(context));
            }
        }
    }

    public Map<String, VariantContext> unmatched() {
        return unmatched;
    }

    public List<StructuralVariant> results() {
        return results;
    }

    public Set<String> sampleNames() {
        return samples;
    }

    private static StructuralVariant create(VariantContext context) {
        final StructuralVariantType type = type(context);
        Preconditions.checkArgument(!StructuralVariantType.BND.equals(type));

        byte startOrientation = 0, endOrientation = 0;
        switch (type) {
            case INV:
                if (context.hasAttribute("INV3")) {
                    startOrientation = 1;
                    endOrientation = -1;
                } else if (context.hasAttribute("INV5")) {
                    startOrientation = -1;
                    endOrientation = 1;
                }
                break;
            case DEL:
                startOrientation = endOrientation = 1;
                break;
            case INS:
                startOrientation = endOrientation = 1;
                break;
            case DUP:
                startOrientation = endOrientation = -1;
                break;
        }

        return ImmutableStructuralVariant.builder()
                .startChromosome(context.getContig())
                .startPosition(context.getStart())
                .startOrientation(startOrientation)
                .startHomology(context.getAttributeAsString(HOM_SEQ, ""))
                .endChromosome(context.getContig())
                .endPosition(context.getEnd())
                .endOrientation(endOrientation)
                .endHomology("")
                .insertSequence(context.getAttributeAsString(INS_SEQ, ""))
                .type(type)
                .build();
    }

    private static StructuralVariant create(VariantContext first, VariantContext second) {
        Preconditions.checkArgument(StructuralVariantType.BND.equals(type(first)));
        Preconditions.checkArgument(StructuralVariantType.BND.equals(type(second)));

        byte startOrientation = 0, endOrientation = 0;
        final String alt = first.getAlternateAllele(0).getDisplayString();
        final String[] leftSplit = alt.split("\\]");
        final String[] rightSplit = alt.split("\\[");
        if (leftSplit.length >= 2) {
            if (leftSplit[0].length() > 0) {
                startOrientation = 1;
                endOrientation = -1;
            } else {
                startOrientation = endOrientation = -1;
            }
        } else if (rightSplit.length >= 2) {
            if (rightSplit[0].length() > 0) {
                startOrientation = endOrientation = 1;
            } else {
                startOrientation = -1;
                endOrientation = 1;
            }
        }

        return ImmutableStructuralVariant.builder()
                .startChromosome(first.getContig())
                .startPosition(first.getStart())
                .startOrientation(startOrientation)
                .startHomology(first.getAttributeAsString(HOM_SEQ, ""))
                .endChromosome(second.getContig())
                .endPosition(second.getEnd())
                .endOrientation(endOrientation)
                .endHomology(second.getAttributeAsString(HOM_SEQ, ""))
                .insertSequence(first.getAttributeAsString(INS_SEQ, ""))
                .type(StructuralVariantType.BND)
                .build();
    }

    private static StructuralVariantType type(VariantContext context) {
        return StructuralVariantType.fromAttribute((String) context.getAttribute(TYPE));
    }
}
