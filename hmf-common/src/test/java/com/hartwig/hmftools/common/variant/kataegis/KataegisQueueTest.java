package com.hartwig.hmftools.common.variant.kataegis;

import static org.junit.Assert.assertEquals;

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
        final VariantContext context1 = create("1", 100, true);
        final VariantContext context2 = create("1", 102, true);
        final VariantContext context3 = create("1", 103, true);
        final VariantContext context4 = create("1", 104, true);
        final VariantContext context5 = create("1", 105, true);
        final VariantContext context6 = create("1", 106, true);

        List<VariantContext> result = kataegis(Lists.newArrayList(context1, context2, context3, context4, context5, context6));
        assertEquals(6, result.size());
        assertEquals(KataegisStatus.FWD.toString(), result.get(0).getAttribute(KataegisQueue.KATAEGIS_FLAG));
    }

    @Test
    public void testChangeContig() {
        final VariantContext context1 = create("1", 100, true);
        final VariantContext context2 = create("1", 102, true);
        final VariantContext context3 = create("1", 103, true);
        final VariantContext context4 = create("2", 104, true);
        final VariantContext context5 = create("2", 105, true);
        final VariantContext context6 = create("2", 106, true);

        List<VariantContext> result = kataegis(Lists.newArrayList(context1, context2, context3, context4, context5, context6));
        assertEquals(6, result.size());
        assertEquals(KataegisStatus.NONE.toString(), result.get(0).getAttribute(KataegisQueue.KATAEGIS_FLAG));
    }

    @Test
    public void testAverageDistanceClustered() {
        final VariantContext context1 = create("1", 100, true);
        final VariantContext context2 = create("1", 102, true);
        final VariantContext context3 = create("1", 103, true);
        final VariantContext context4 = create("1", 104, true);
        final VariantContext context5 = create("1", 105, true);
        final VariantContext context6 = create("1", 5100, true);

        List<VariantContext> result = kataegis(Lists.newArrayList(context1, context2, context3, context4, context5, context6));
        assertEquals(6, result.size());
        assertEquals(KataegisStatus.FWD.toString(), result.get(0).getAttribute(KataegisQueue.KATAEGIS_FLAG));
    }

    @Test
    public void testAverageDistanceSpread() {
        final VariantContext context1 = create("1", 1100, true);
        final VariantContext context2 = create("1", 2102, true);
        final VariantContext context3 = create("1", 3103, true);
        final VariantContext context4 = create("1", 4104, true);
        final VariantContext context5 = create("1", 5105, true);
        final VariantContext context6 = create("1", 6104, true);

        List<VariantContext> result = kataegis(Lists.newArrayList(context1, context2, context3, context4, context5, context6));
        assertEquals(6, result.size());
        assertEquals(KataegisStatus.FWD.toString(), result.get(0).getAttribute(KataegisQueue.KATAEGIS_FLAG));
    }

    @Test
    public void testMaxDepthAcceptable() {
        final VariantContext context1 = create("1", 100, true);
        final VariantContext context2 = create("1", 102, true);
        final VariantContext context3 = create("1", 103, true);
        final VariantContext context4 = create("1", 104, true);
        final VariantContext context5 = create("1", 105, true);
        final VariantContext context6 = create("1", 5106, true);

        List<VariantContext> result = kataegis(Lists.newArrayList(context1, context2, context3, context4, context5, context6));
        assertEquals(6, result.size());
        assertEquals(KataegisStatus.NONE.toString(), result.get(0).getAttribute(KataegisQueue.KATAEGIS_FLAG));
    }

    @NotNull
    static VariantContext create(String contig, long start, boolean kataegis) {
        return create(contig, start, "C", kataegis ? "T" : "A");
    }

    @NotNull
    private static VariantContext create(String contig, long start, @NotNull final String ref, @NotNull final String alt) {
        Allele refAllele = Allele.create(ref, true);
        Allele altAllele = Allele.create(alt, false);

        return new VariantContextBuilder("Source", contig, start, start, Lists.newArrayList(refAllele, altAllele)).make();
    }

    @NotNull
    private static List<VariantContext> kataegis(@NotNull final List<VariantContext> contexts) {
        final Predicate<VariantContext> kataegisPredicate = context -> context.getAlternateAllele(0).getBaseString().equals("T");
        final List<VariantContext> result = Lists.newArrayList();
        KataegisQueue inner = new KataegisQueue(KataegisStatus.FWD, kataegisPredicate, result::add);
        contexts.forEach(inner::accept);
        inner.flush();
        return result;
    }
}
