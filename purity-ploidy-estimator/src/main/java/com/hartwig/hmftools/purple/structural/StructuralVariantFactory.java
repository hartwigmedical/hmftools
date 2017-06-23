package com.hartwig.hmftools.purple.structural;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructualVariantType;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

class StructuralVariantFactory {

    private static String TYPE = "SVTYPE";
    private static String MATE_ID = "MATEID";

    private final Map<String, VariantContext> unmatched = new HashMap<>();
    private final List<StructuralVariant> results = Lists.newArrayList();
    private final PassingVariantFilter filter = new PassingVariantFilter();
    private final Set<String> samples = Sets.newHashSet();

    void addVariantContext(VariantContext context) {
        if (filter.test(context)) {

            samples.addAll(context.getSampleNames());
            final StructualVariantType type = type(context);
            if (type.equals(StructualVariantType.BND)) {
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
        final StructualVariantType type = type(context);
        Preconditions.checkArgument(!StructualVariantType.BND.equals(type));

        return ImmutableStructuralVariant.builder()
                .startChromosome(context.getContig())
                .startPosition(context.getStart())
                .endChromosome(context.getContig())
                .endPosition(context.getEnd())
                .type(type)
                .build();
    }

    private static StructuralVariant create(VariantContext first, VariantContext second) {
        Preconditions.checkArgument(StructualVariantType.BND.equals(type(first)));
        Preconditions.checkArgument(StructualVariantType.BND.equals(type(second)));

        return ImmutableStructuralVariant.builder()
                .startChromosome(first.getContig())
                .startPosition(first.getStart())
                .endChromosome(second.getContig())
                .endPosition(second.getEnd())
                .type(StructualVariantType.BND)
                .build();
    }

    private static StructualVariantType type(VariantContext context) {
        return StructualVariantType.fromAttribute((String) context.getAttribute(TYPE));
    }
}
