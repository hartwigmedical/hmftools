package com.hartwig.hmftools.sage;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class SageHotspotAnnotationTest {

    private static final String SAMPLE = "sample";

    private VCFCodec codec;

    @Before
    public void setup() {
        codec = createTestCodec();
    }

    @Test
    public void testComparatorUsesAlleles() {
        final String line1 = "15\t12345678\trs1;UCSC\tC\tA,G\t2\tPASS\tinfo;\tGT:AD:DP\t0/1:60,60:121";
        final String line2 = "15\t12345678\trs1;UCSC\tC\tA,GT\t2\tPASS\tinfo;\tGT:AD:DP\t0/1:60,60:121";

        SAMSequenceDictionary samSequenceDictionary = new SAMSequenceDictionary();
        samSequenceDictionary.addSequence(new SAMSequenceRecord("15", 123456789));

        SageHotspotAnnotation.VCComparator comparator = new SageHotspotAnnotation.VCComparator(samSequenceDictionary);
        assertEquals(0, comparator.compare(codec.decode(line1), codec.decode(line1)));
        assertNotEquals(0, comparator.compare(codec.decode(line1), codec.decode(line2)));
    }

    @NotNull
    private static VCFCodec createTestCodec() {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet(SAMPLE));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }
}
