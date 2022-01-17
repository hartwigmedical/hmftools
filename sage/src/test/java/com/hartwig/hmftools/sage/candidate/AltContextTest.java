package com.hartwig.hmftools.sage.candidate;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.ReadContextTest;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class AltContextTest
{

    private static final String CHROM = "1";
    private static final int POS = 1000;

    @Test
    public void testIncompleteReadContext()
    {
        final RefContext refContext = new RefContext(CHROM, POS, false);
        final AltContext victim = new AltContext(refContext, "C", "T");
        final String core1 = "GATAC";
        final String core2 = "GATAA";

        final List<ReadContext> readContexts = Lists.newArrayList();
        readContexts.add(simpleSnv("A", core1, "AG"));
        readContexts.add(simpleSnv("A", core1, "AG"));
        readContexts.add(simpleSnv("T", core2, "TT"));

        // Ordering should not matter!
        Collections.shuffle(readContexts);
        readContexts.forEach(x -> victim.addReadContext(0, x));

        assertFalse(victim.finaliseAndValidate());
    }

    @Test
    public void testFullMatch()
    {
        final RefContext refContext = new RefContext(CHROM, POS, false);
        final AltContext victim = new AltContext(refContext, "C", "T");
        final String core1 = "GATAC";

        final List<ReadContext> readContexts = Lists.newArrayList();
        readContexts.add(simpleSnv("AG", core1, "AG"));
        readContexts.add(simpleSnv("AG", core1, "AG"));
        readContexts.add(simpleSnv("TT", core1, "TT"));

        // Ordering should not matter!
        Collections.shuffle(readContexts);
        readContexts.forEach(x -> victim.addReadContext(0, x));

        assertTrue(victim.finaliseAndValidate());
        assertEquals("AG" + core1 + "AG", new String(victim.readContext().readBases()));
    }

    @Test
    public void testCoreMatchAfterFullMatch()
    {
        final RefContext refContext = new RefContext(CHROM, POS, false);
        final AltContext victim = new AltContext(refContext, "C", "T");

        final String core1 = "GATAC";
        final String core2 = "GATAA";
        assertNotEquals(core1, core2);

        final List<ReadContext> readContexts = Lists.newArrayList();
        // At this point either is valid
        readContexts.add(simpleSnv("AG", core2, "AG"));
        readContexts.add(simpleSnv("AG", core2, "AG"));

        readContexts.add(simpleSnv("AG", core1, "AG"));
        readContexts.add(simpleSnv("AG", core1, "AG"));

        // Add new core1 (with different flanks)
        readContexts.add(simpleSnv("TT", core1, "TT"));

        // Ordering should not matter!
        Collections.shuffle(readContexts);
        readContexts.forEach(x -> victim.addReadContext(0, x));

        assertTrue(victim.finaliseAndValidate());
        assertEquals("AG" + core1 + "AG", new String(victim.readContext().readBases()));
    }

    @Test
    public void testPartialMatchAfterFullMatch()
    {
        final RefContext refContext = new RefContext(CHROM, POS, false);
        final AltContext victim = new AltContext(refContext, "C", "T");

        final String core1 = "GATAC";
        final String core2 = "GATAA";
        assertNotEquals(core1, core2);

        final List<ReadContext> readContexts = Lists.newArrayList();
        // At this point either is valid
        readContexts.add(simpleSnv("AG", core2, "AG"));
        readContexts.add(simpleSnv("AG", core2, "AG"));

        readContexts.add(simpleSnv("AG", core1, "AG"));
        readContexts.add(simpleSnv("AG", core1, "AG"));

        // Add new core1 (with different flanks)
        readContexts.add(simpleSnv("", core1, "AG"));

        // Ordering should not matter!
        Collections.shuffle(readContexts);
        readContexts.forEach(x -> victim.addReadContext(0, x));

        assertTrue(victim.finaliseAndValidate());
        assertEquals("AG" + core1 + "AG", new String(victim.readContext().readBases()));
    }

    @Test
    public void testBalancedCore()
    {
        String core = "CAT";

        final RefContext refContext = new RefContext(CHROM, POS, false);
        final AltContext victim = new AltContext(refContext, "C", "T");

        victim.addReadContext(3, simpleReadContext("AAA", core, "CCC"));
        assertEquals("AAACATCCC", new String(victim.interimReadContexts().get(0).readContext().readBases()));

        // Adding same size shouldn't change it
        victim.addReadContext(2, simpleReadContext("TTA", core, "CGG"));
        assertEquals("AAACATCCC", new String(victim.interimReadContexts().get(0).readContext().readBases()));

        // Adding one to left won't change it
        victim.addReadContext(3, simpleReadContext("TAAA", core, "CCC"));
        assertEquals("AAACATCCC", new String(victim.interimReadContexts().get(0).readContext().readBases()));

        // Adding one to right won't change it
        victim.addReadContext(3, simpleReadContext("AAA", core, "CCCG"));
        assertEquals("AAACATCCC", new String(victim.interimReadContexts().get(0).readContext().readBases()));

        // Adding one to both WILL change it
        victim.addReadContext(3, simpleReadContext("TAAA", core, "CCCG"));
        AltContext.ReadContextCandidate finalCandidate = victim.interimReadContexts().get(0);

        assertEquals("TAAACATCCCG", new String(finalCandidate.readContext().readBases()));
        assertEquals(2, finalCandidate.minNumberOfEvents());
        assertEquals(4, finalCandidate.fullMatch());
    }

    @NotNull
    public static ReadContext simpleReadContext(@NotNull final String leftFlank, @NotNull final String core,
            @NotNull final String rightFlank)
    {
        assert (core.length() == 3);
        return ReadContextTest.create(POS, 1, leftFlank, core, rightFlank, 1);
    }

    @NotNull
    public static ReadContext simpleSnv(@NotNull final String leftFlank, @NotNull final String core, @NotNull final String rightFlank)
    {
        assert (core.length() == 5);
        return ReadContextTest.simpleSnv(POS, leftFlank, core, rightFlank);
    }
}
