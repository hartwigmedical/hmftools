package com.hartwig.hmftools.amber.contamination;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.amber.VafReading;
import com.hartwig.hmftools.common.amber.AmberBase;
import com.hartwig.hmftools.common.amber.BaseDepthData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Test;

public class BinomialVafPredicateTest
{
    @Test
    public void testTest()
    {
        // Sigma = sqrt(n*p*(1-p))
        // 100*0.2*0.8 = 16 -> 4 -> 8 -> [12,28]
        BinomialVafPredicate pred = new BinomialVafPredicate(0.2);
        assertFalse(pred.test(tc(100, 10, 5)));
        assertFalse(pred.test(tc(100, 10, 10)));
        assertFalse(pred.test(tc(100, 10, 11)));
        assertTrue(pred.test(tc(100, 10, 12)));
        assertTrue(pred.test(tc(100, 10, 13)));
        assertTrue(pred.test(tc(100, 10, 21)));
        assertTrue(pred.test(tc(100, 10, 25)));
        assertTrue(pred.test(tc(100, 10, 26)));
        assertTrue(pred.test(tc(100, 10, 27)));
        assertTrue(pred.test(tc(100, 10, 28)));
        assertFalse(pred.test(tc(100, 10, 29)));
        assertFalse(pred.test(tc(100, 10, 30)));

        // 1000*.1 *.9 = 90 -> 9.487 -> 18.97 -> [82,118] (as integers)
        pred = new BinomialVafPredicate(0.1);
        assertFalse(pred.test(tc(1000, 900, 5)));
    }

    private VafReading tc(int readDepth, int refSupport, int altSupport)
    {
        return new VafReading(V38.versionedChromosome(HumanChromosome._2), 100_000_000, readDepth, refSupport, altSupport);
    }
}
