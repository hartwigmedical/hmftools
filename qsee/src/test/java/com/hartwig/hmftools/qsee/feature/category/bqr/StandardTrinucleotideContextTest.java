package com.hartwig.hmftools.qsee.feature.category.bqr;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class StandardTrinucleotideContextTest
{
    @Test
    public void canStandardiseTtoGContext()
    {
        StandardTrinucleotideContext context;

        context = StandardTrinucleotideContext.fromNucleotides('T', 'G', "GTA");
        assertEquals("G[T>G]A +", context.mutationContextString());

        context = StandardTrinucleotideContext.fromNucleotides('A', 'C', "TAC");
        assertEquals("G[T>G]A -", context.mutationContextString());
    }

    @Test
    public void canStandardiseCtoTContext()
    {
        StandardTrinucleotideContext context;

        context = StandardTrinucleotideContext.fromNucleotides('C', 'T', "GCA");
        assertEquals("G[C>T]A +", context.mutationContextString());

        context = StandardTrinucleotideContext.fromNucleotides('G', 'A', "TGC");
        assertEquals("G[C>T]A -", context.mutationContextString());
    }
}
