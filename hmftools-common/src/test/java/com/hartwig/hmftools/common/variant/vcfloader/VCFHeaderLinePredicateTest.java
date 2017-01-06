package com.hartwig.hmftools.common.variant.vcfloader;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class VCFHeaderLinePredicateTest {

    @Test
    public void vcfHeaderLinePredicateWorks() {
        VCFHeaderLinePredicate predicate = new VCFHeaderLinePredicate();
        assertTrue(predicate.test("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample"));
        assertFalse(predicate.test("##fileformat=VCFv4.1"));
        assertFalse(predicate.test("1\t1\t.\tC\tA\t.\tPASS\tset=freebayes\tGT:AD:DP\t0/1:48,15:64"));
    }
}