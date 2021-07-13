package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.isofox.novel.cohort.AcceptorDonorType.ACCEPTOR;
import static com.hartwig.hmftools.isofox.novel.cohort.AcceptorDonorType.DONOR;

import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.isofox.novel.cohort.SpliceVariant;

import org.junit.Test;

public class SpliceVariantMatchingTest
{
    @Test
    public void testSpliceSiteBaseContext()
    {
        final String chromosome = "1";

        final MockRefGenome refGenome = new MockRefGenome();

        //                       01234567890123456789012345678901234567890
        final String refBases = "AGCTAGCTAGAGCTAGCTAGCTAGCTAGCTAAGCTAGCTAG";
        refGenome.RefGenomeMap.put(chromosome, refBases);

        // SNV at the same base
        String baseContext = SpliceVariant.getBaseContext(
                chromosome, 20, "C", "G", 20, DONOR, refGenome);

        assertTrue(baseContext.equals("GGTAGCTAGCTA"));

        // SNV at an earlier base
        baseContext = SpliceVariant.getBaseContext(
                chromosome, 18, "A", "G", 20, ACCEPTOR, refGenome);

        assertTrue(baseContext.equals("AGCTAGCTGGCT"));

        // SNV at a later base
        baseContext = SpliceVariant.getBaseContext(
                chromosome, 25, "T", "G", 20, DONOR, refGenome);

        assertTrue(baseContext.equals("GCTAGCGAGCTA"));

        // DEL at the same base
        baseContext = SpliceVariant.getBaseContext(
                chromosome, 20, "CTA", "C", 20, DONOR, refGenome);

        assertTrue(baseContext.equals("GCGCTAGCTAAG"));

        // DEL at a later base
        baseContext = SpliceVariant.getBaseContext(
                chromosome, 22, "AGC", "A", 20, DONOR, refGenome);

        assertTrue(baseContext.equals("GCTATAGCTAAG"));

        // DEL at an earlier base
        baseContext = SpliceVariant.getBaseContext(
                chromosome, 15, "GCT", "G", 20, ACCEPTOR, refGenome);

        assertTrue(baseContext.equals("AGAGCTAGAGCT"));

        // INS at the same base
        baseContext = SpliceVariant.getBaseContext(
                chromosome, 20, "C", "CCCCC", 20, DONOR, refGenome);

        assertTrue(baseContext.equals("GCCCCCTAGCTA"));

        // INS at a later base
        baseContext = SpliceVariant.getBaseContext(
                chromosome, 25, "T", "TTTTT", 20, DONOR, refGenome);

        assertTrue(baseContext.equals("GCTAGCTTTTTA"));

        // INS at same base longer than required
        baseContext = SpliceVariant.getBaseContext(
                chromosome, 20, "C", "CCCCC", 20, ACCEPTOR, refGenome);

        assertTrue(baseContext.equals("AGCTAGCTAGCC"));

        // INS at a later base longer than required
        baseContext = SpliceVariant.getBaseContext(
                chromosome, 21, "T", "TTT", 20, ACCEPTOR, refGenome);

        assertTrue(baseContext.equals("AGCTAGCTAGCT"));

        // INS at an earlier base
        baseContext = SpliceVariant.getBaseContext(
                chromosome, 18, "A", "AAAAA", 20, ACCEPTOR, refGenome);

        assertTrue(baseContext.equals("AGCTAAAAAGCT"));

        // variants outside the required context range have no effect on the base context
        baseContext = SpliceVariant.getBaseContext(
                chromosome, 23, "A", "G", 20, ACCEPTOR, refGenome);

        assertTrue(baseContext.equals("AGCTAGCTAGCT"));

        baseContext = SpliceVariant.getBaseContext(
                chromosome, 17, "A", "G", 20, DONOR, refGenome);

        assertTrue(baseContext.equals("GCTAGCTAGCTA"));
    }
}
