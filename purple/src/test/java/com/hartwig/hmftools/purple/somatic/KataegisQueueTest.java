package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.KATAEGIS_FLAG;
import static com.hartwig.hmftools.purple.MiscTestUtils.SAMPLE_ID;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Predicate;

import com.google.common.collect.Lists;

import org.junit.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class KataegisQueueTest
{
    @Test
    public void testExpectedBehaviour()
    {
        final List<SomaticVariant> input = Lists.newArrayList();
        for(int i = 0; i < KataegisQueue.MIN_COUNT; i++)
        {
            input.add(create("1", 100 + i, true));
        }

        List<SomaticVariant> result = kataegis(input);
        assertEquals(KataegisQueue.MIN_COUNT, result.size());
        assertEquals("TST_2", result.get(0).context().getAttribute(KATAEGIS_FLAG));
    }

    @Test
    public void testChangeContig()
    {
        final List<SomaticVariant> input = Lists.newArrayList();
        for(int i = 0; i < KataegisQueue.MIN_COUNT; i++)
        {
            input.add(create(i % 2 == 0 ? "1" : "2", 100 + i, true));
        }

        List<SomaticVariant> result = kataegis(input);
        assertEquals(KataegisQueue.MIN_COUNT, result.size());
        assertFalse(result.get(0).context().hasAttribute(KATAEGIS_FLAG));
    }

    @Test
    public void testAverageDistanceClustered()
    {
        final SomaticVariant context1 = create("1", 100, true);
        final SomaticVariant context2 = create("1", 102, true);
        final SomaticVariant context3 = create("1", 103, true);
        final SomaticVariant context4 = create("1", 104, true);
        final SomaticVariant context5 = create("1", 105, true);
        final SomaticVariant context6 = create("1", 2100, true);

        List<SomaticVariant> result = kataegis(Lists.newArrayList(context1, context2, context3, context4, context5, context6));
        assertEquals(6, result.size());
        assertEquals("TST_2", result.get(0).context().getAttribute(KATAEGIS_FLAG));
    }

    @Test
    public void testAverageDistanceSpread()
    {
        final SomaticVariant context1 = create("1", 1101, true);
        final SomaticVariant context2 = create("1", 2102, true);
        final SomaticVariant context3 = create("1", 3103, true);
        final SomaticVariant context4 = create("1", 4104, true);
        final SomaticVariant context5 = create("1", 5105, true);
        final SomaticVariant context6 = create("1", 6103, true);

        List<SomaticVariant> result = kataegis(Lists.newArrayList(context1, context2, context3, context4, context5, context6));
        assertEquals(6, result.size());
        assertEquals("TST_2", result.get(0).context().getAttribute(KATAEGIS_FLAG));
    }

    @Test
    public void testMaxDepthAcceptable()
    {
        final List<SomaticVariant> input = Lists.newArrayList();
        int i;
        for(i = 0; i < KataegisQueue.MIN_COUNT; i++)
        {
            input.add(create("1", 100 + i, true));
        }
        i--;

        input.add(create("1", 100 + i + KataegisQueue.MAX_ABS_DISTANCE, true));
        input.add(create("1", 100 + i + 2 * KataegisQueue.MAX_ABS_DISTANCE + 1, true));

        final List<SomaticVariant> result = kataegis(input);
        assertEquals(KataegisQueue.MIN_COUNT + 2, result.size());
        assertEquals("TST_2", result.get(i + 1).context().getAttribute(KATAEGIS_FLAG));
        assertFalse(result.get(i + 2).context().hasAttribute(KATAEGIS_FLAG));
    }

    static SomaticVariant create(final String contig, long start, boolean kataegis)
    {
        return create(contig, start, kataegis ? "T" : "A");
    }

    private static SomaticVariant create(final String contig, long start, final String alt)
    {
        Allele refAllele = Allele.create("C", true);
        Allele altAllele = Allele.create(alt, false);

        VariantContext context = new VariantContextBuilder("Source", contig, start, start, Lists.newArrayList(refAllele, altAllele)).make();
        return new SomaticVariant(context, SAMPLE_ID, null);
    }

    private static List<SomaticVariant> kataegis(final List<SomaticVariant> variants)
    {
        final Predicate<SomaticVariant> kataegisPredicate = variant -> variant.context().getAlternateAllele(0).getBaseString().equals("T");
        final List<SomaticVariant> result = Lists.newArrayList();
        KataegisQueue inner = new KataegisQueue("TST", new AtomicInteger(), kataegisPredicate, result::add);
        variants.forEach(inner::processVariant);
        inner.flush();
        return result;
    }
}
