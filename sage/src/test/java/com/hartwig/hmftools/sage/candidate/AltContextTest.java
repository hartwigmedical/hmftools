package com.hartwig.hmftools.sage.candidate;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;

import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.junit.Test;

public class AltContextTest
{
    private static final int POSITION = 1000;
    private static final String REF = "C";
    private static final String ALT = "T";

    private static VariantReadContext simpleSnv(final String leftFlank, final String core, final String rightFlank)
    {
        SimpleVariant variant = createSimpleVariant(POSITION, REF, ALT);
        String leftCore = core.substring(0, 2);
        String rightCore = core.substring(3, 5);
        return createReadContext(variant, leftCore, rightCore, leftFlank, rightFlank);
    }

    /* TODO: convert to read-based tests
    @Test
    public void testIncompleteReadContext()
    {
        RefContext refContext = new RefContext(CHR_1, POSITION);
        AltContext altContext = new AltContext(refContext, "C", "T");
        final String core1 = "GATAC";
        final String core2 = "GATAA";

        final List<VariantReadContext> readContexts = Lists.newArrayList();

        readContexts.add(simpleSnv("A", core1, "AG"));
        readContexts.add(simpleSnv("A", core1, "AG"));
        readContexts.add(simpleSnv("T", core2, "TT"));

        readContexts.forEach(x -> altContext.addReadContext(0, x));

        // this test is a misnomer - the above read contexts are all valid since there is no incomplete flank concept anymore
        altContext.selectCandidates();
        assertTrue(altContext.hasValidCandidate());
    }

    @Test
    public void testFullMatch()
    {
        final RefContext refContext = new RefContext(CHR_1, POSITION);
        final AltContext altContext = new AltContext(refContext, "C", "T");
        final String core1 = "GATAC";

        final List<VariantReadContext> readContexts = Lists.newArrayList();
        readContexts.add(simpleSnv("AG", core1, "AG"));
        readContexts.add(simpleSnv("AG", core1, "AG"));
        readContexts.add(simpleSnv("TT", core1, "TT"));

        // Ordering should not matter!
        Collections.shuffle(readContexts);
        readContexts.forEach(x -> altContext.addReadContext(0, x));

        altContext.selectCandidates();
        assertTrue(altContext.hasValidCandidate());
        assertEquals("AG" + core1 + "AG", new String(altContext.readContext().readBases()));
    }

    @Test
    public void testCoreMatchAfterFullMatch()
    {
        final RefContext refContext = new RefContext(CHR_1, POSITION);
        final AltContext altContext = new AltContext(refContext, "C", "T");

        final String core1 = "GATAC";
        final String core2 = "GATAA";

        final List<VariantReadContext> readContexts = Lists.newArrayList();
        // At this point either is valid
        readContexts.add(simpleSnv("AG", core2, "AG"));
        readContexts.add(simpleSnv("AG", core2, "AG"));

        readContexts.add(simpleSnv("AG", core1, "AG"));
        readContexts.add(simpleSnv("AG", core1, "AG"));

        // Add new core1 (with different flanks)
        readContexts.add(simpleSnv("TT", core1, "TT"));

        readContexts.forEach(x -> altContext.addReadContext(0, x));

        altContext.selectCandidates();
        assertTrue(altContext.hasValidCandidate());
        assertEquals("AG" + core1 + "AG", new String(altContext.readContext().readBases()));
    }
    */
}
