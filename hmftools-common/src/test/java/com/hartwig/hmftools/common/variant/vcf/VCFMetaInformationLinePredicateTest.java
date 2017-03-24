package com.hartwig.hmftools.common.variant.vcf;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class VCFMetaInformationLinePredicateTest {

    @Test
    public void vcfMetaInformationLinePredicateWorks() {
        final VCFMetaInformationLinePredicate predicate = new VCFMetaInformationLinePredicate();
        assertFalse(predicate.test("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample"));
        assertTrue(predicate.test("##fileformat=VCFv4.1"));
        assertFalse(predicate.test("1\t1\t.\tC\tA\t.\tPASS\tset=freebayes\tGT:AD:DP\t0/1:48,15:64"));
    }
}