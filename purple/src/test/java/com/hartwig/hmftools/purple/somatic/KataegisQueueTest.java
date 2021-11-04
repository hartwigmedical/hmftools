package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.KATAEGIS_FLAG;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.List;
import java.util.function.Predicate;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class KataegisQueueTest {

    @Test
    public void testExpectedBehaviour() {
        final List<VariantContext> input = Lists.newArrayList();
        for (int i = 0; i < KataegisQueue.MIN_COUNT; i++) {
            input.add(create("1", 100 + i, true));
        }

        List<VariantContext> result = kataegis(input);
        assertEquals(KataegisQueue.MIN_COUNT, result.size());
        assertEquals("TST_1", result.get(0).getAttribute(KATAEGIS_FLAG));
    }

    @Test
    public void testChangeContig() {
        final List<VariantContext> input = Lists.newArrayList();
        for (int i = 0; i < KataegisQueue.MIN_COUNT; i++) {
            input.add(create(i % 2 == 0 ? "1" : "2", 100 + i, true));
        }

        List<VariantContext> result = kataegis(input);
        assertEquals(KataegisQueue.MIN_COUNT, result.size());
        assertFalse(result.get(0).hasAttribute(KATAEGIS_FLAG));
    }

    @Test
    public void testAverageDistanceClustered() {
        final VariantContext context1 = create("1", 100, true);
        final VariantContext context2 = create("1", 102, true);
        final VariantContext context3 = create("1", 103, true);
        final VariantContext context4 = create("1", 104, true);
        final VariantContext context5 = create("1", 105, true);
        final VariantContext context6 = create("1", 2100, true);

        List<VariantContext> result = kataegis(Lists.newArrayList(context1, context2, context3, context4, context5, context6));
        assertEquals(6, result.size());
        assertEquals("TST_1", result.get(0).getAttribute(KATAEGIS_FLAG));
    }

    @Test
    public void testAverageDistanceSpread() {
        final VariantContext context1 = create("1", 1101, true);
        final VariantContext context2 = create("1", 2102, true);
        final VariantContext context3 = create("1", 3103, true);
        final VariantContext context4 = create("1", 4104, true);
        final VariantContext context5 = create("1", 5105, true);
        final VariantContext context6 = create("1", 6103, true);

        List<VariantContext> result = kataegis(Lists.newArrayList(context1, context2, context3, context4, context5, context6));
        assertEquals(6, result.size());
        assertEquals("TST_1", result.get(0).getAttribute(KATAEGIS_FLAG));
    }

    @Test
    public void testMaxDepthAcceptable() {
        final List<VariantContext> input = Lists.newArrayList();
        int i;
        for (i = 0; i < KataegisQueue.MIN_COUNT; i++) {
            input.add(create("1", 100 + i, true));
        }
        i--;

        input.add(create("1", 100 + i + KataegisQueue.MAX_ABS_DISTANCE, true));
        input.add(create("1", 100 + i + 2 * KataegisQueue.MAX_ABS_DISTANCE + 1, true));

        final List<VariantContext> result = kataegis(input);
        assertEquals(KataegisQueue.MIN_COUNT + 2, result.size());
        assertEquals("TST_1", result.get(i + 1).getAttribute(KATAEGIS_FLAG));
        assertFalse(result.get(i + 2).hasAttribute(KATAEGIS_FLAG));
    }

    @NotNull
    static VariantContext create(@NotNull String contig, long start, boolean kataegis) {
        return create(contig, start, kataegis ? "T" : "A");
    }

    @NotNull
    private static VariantContext create(@NotNull String contig, long start, @NotNull final String alt) {
        Allele refAllele = Allele.create("C", true);
        Allele altAllele = Allele.create(alt, false);

        return new VariantContextBuilder("Source", contig, start, start, Lists.newArrayList(refAllele, altAllele)).make();
    }

    @NotNull
    private static List<VariantContext> kataegis(@NotNull final List<VariantContext> contexts) {
        final Predicate<VariantContext> kataegisPredicate = context -> context.getAlternateAllele(0).getBaseString().equals("T");
        final List<VariantContext> result = Lists.newArrayList();
        KataegisQueue inner = new KataegisQueue("TST", kataegisPredicate, result::add);
        contexts.forEach(inner::accept);
        inner.flush();
        return result;
    }
}
