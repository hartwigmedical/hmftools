package com.hartwig.hmftools.common.variant.structural;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public class StructuralVariantFactory {

    private final static String TYPE = "SVTYPE";
    private final static String MATE_ID = "MATEID";
    private final static String INS_SEQ = "SVINSSEQ";
    private final static String HOM_SEQ = "HOMSEQ";
    private final static String BPI_START = "BPI_START";
    private final static String BPI_END = "BPI_END";
    private final static String BPI_AF = "BPI_AF";

    private final Map<String, VariantContext> unmatched = new HashMap<>();
    private final List<StructuralVariant> results = Lists.newArrayList();
    private final PassingVariantFilter filter = new PassingVariantFilter();
    private final Set<String> samples = Sets.newHashSet();

    public void addVariantContext(VariantContext context) {
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

        final int start = context.hasAttribute(BPI_START) ? context.getAttributeAsInt(BPI_START, -1) : context.getStart();
        final int end = context.hasAttribute(BPI_END) ? context.getAttributeAsInt(BPI_END, -1) : context.getEnd();
        final List<Double> af = context.hasAttribute(BPI_AF) ? context.getAttributeAsDoubleList(BPI_AF, 0.0) : Collections.emptyList();

        byte startOrientation = 0, endOrientation = 0;
        switch (type) {
            case INV:
                if (context.hasAttribute("INV3")) {
                    startOrientation = endOrientation = 1;
                } else if (context.hasAttribute("INV5")) {
                    startOrientation = endOrientation = -1;
                }
                break;
            case DEL:
                startOrientation = 1;
                endOrientation = -1;
                break;
            case INS:
                startOrientation = 1;
                endOrientation = -1;
                break;
            case DUP:
                startOrientation = -1;
                endOrientation = 1;
                break;
        }

        return ImmutableStructuralVariant.builder()
                .startChromosome(context.getContig())
                .startPosition(start)
                .startOrientation(startOrientation)
                .startHomology(context.getAttributeAsString(HOM_SEQ, ""))
                .startAF(af.size() == 2 ? af.get(0) : null)
                .endChromosome(context.getContig())
                .endPosition(end)
                .endOrientation(endOrientation)
                .endHomology("")
                .endAF(af.size() == 2 ? af.get(1) : null)
                .insertSequence(context.getAttributeAsString(INS_SEQ, ""))
                .type(type)
                .build();
    }

    private static StructuralVariant create(VariantContext first, VariantContext second) {
        Preconditions.checkArgument(StructuralVariantType.BND.equals(type(first)));
        Preconditions.checkArgument(StructuralVariantType.BND.equals(type(second)));

        final int start = first.hasAttribute(BPI_START) ? first.getAttributeAsInt(BPI_START, -1) : first.getStart();
        final int end = second.hasAttribute(BPI_START) ? second.getAttributeAsInt(BPI_START, -1) : second.getStart();
        final List<Double> af = first.hasAttribute(BPI_AF) ? first.getAttributeAsDoubleList(BPI_AF, 0.0) : Collections.emptyList();

        byte startOrientation = 0, endOrientation = 0;
        final String alt = first.getAlternateAllele(0).getDisplayString();
        final String[] leftSplit = alt.split("\\]");
        final String[] rightSplit = alt.split("\\[");
        if (leftSplit.length >= 2) {
            if (leftSplit[0].length() > 0) {
                startOrientation = endOrientation = 1;
            } else {
                startOrientation = -1;
                endOrientation = 1;
            }
        } else if (rightSplit.length >= 2) {
            if (rightSplit[0].length() > 0) {
                startOrientation = 1;
                endOrientation = -1;
            } else {
                startOrientation = endOrientation = -1;
            }
        }

        return ImmutableStructuralVariant.builder()
                .startChromosome(first.getContig())
                .startPosition(start)
                .startOrientation(startOrientation)
                .startHomology(first.getAttributeAsString(HOM_SEQ, ""))
                .startAF(af.size() == 2 ? af.get(0) : null)
                .endChromosome(second.getContig())
                .endPosition(end)
                .endOrientation(endOrientation)
                .endHomology(second.getAttributeAsString(HOM_SEQ, ""))
                .endAF(af.size() == 2 ? af.get(1) : null)
                .insertSequence(first.getAttributeAsString(INS_SEQ, ""))
                .type(StructuralVariantType.BND)
                .build();
    }

    private static StructuralVariantType type(VariantContext context) {
        return StructuralVariantType.fromAttribute((String) context.getAttribute(TYPE));
    }
}
