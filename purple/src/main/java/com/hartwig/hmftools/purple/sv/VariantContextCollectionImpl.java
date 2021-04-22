package com.hartwig.hmftools.purple.sv;

import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;
import java.util.function.Predicate;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFHeader;

public class VariantContextCollectionImpl implements VariantContextCollection {

    private final TreeSet<VariantContext> variantContexts;
    private final List<StructuralVariant> variants = Lists.newArrayList();

    private boolean modified;

    public VariantContextCollectionImpl(@NotNull final VCFHeader header) {
        variantContexts = new TreeSet<>(new VCComparator(header.getSequenceDictionary()));
    }

    @VisibleForTesting
    VariantContextCollectionImpl(@NotNull final List<String> contigs) {
        variantContexts = new TreeSet<>(new VCComparator(contigs));
    }

    @Override
    public void add(@NotNull final VariantContext variantContext) {
        modified = true;
        if (variantContext.contains(variantContext)) {
            variantContexts.remove(variantContext);
            variantContexts.add(variantContext);
        }
    }

    @Override
    public int remove(@NotNull final Predicate<VariantContext> removePredicate) {
        int removed = 0;
        final Iterator<VariantContext> iterator = variantContexts.iterator();
        while (iterator.hasNext()) {
            final VariantContext variantContext = iterator.next();
            if (removePredicate.test(variantContext)) {
                iterator.remove();
                removed++;
                modified = true;
            }
        }
        return removed;
    }

    @NotNull
    @Override
    public List<StructuralVariant> segmentationVariants() {
        if (modified) {
            modified = false;
            final StructuralVariantFactory factory = new StructuralVariantFactory(new SegmentationVariantsFilter());
            variantContexts.forEach(factory::addVariantContext);

            variants.clear();
            variants.addAll(factory.results());
        }

        return variants;
    }

    @NotNull
    @Override
    public Iterator<VariantContext> iterator() {
        return variantContexts.iterator();
    }

    private static class VCComparator extends VariantContextComparator {

        VCComparator(final List<String> contigs) {
            super(contigs);
        }

        VCComparator(final SAMSequenceDictionary dictionary) {
            super(dictionary);
        }

        @Override
        public int compare(final VariantContext firstVariantContext, final VariantContext secondVariantContext) {
            int positionResult = super.compare(firstVariantContext, secondVariantContext);

            return positionResult == 0 ? firstVariantContext.getID().compareTo(secondVariantContext.getID()) : positionResult;
        }
    }
}
