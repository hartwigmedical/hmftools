package com.hartwig.hmftools.serve.sources.vicc.extractors;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.sources.vicc.extractor.GeneRangeExtractor;

import org.junit.Test;

public class GeneRangeExtractorTest {

    @Test
    public void canExtractCodon() {
        assertEquals(Integer.valueOf("600"), GeneRangeExtractor.extractCodonNumber("BRAF (V600)"));
        assertEquals(Integer.valueOf("742"), GeneRangeExtractor.extractCodonNumber("W742"));
        assertEquals(Integer.valueOf("179"), GeneRangeExtractor.extractCodonNumber("Q179X"));
        assertEquals(Integer.valueOf("61"), GeneRangeExtractor.extractCodonNumber("KRAS Q61X"));
    }

    @Test
    public void canExtractExon() {
        assertEquals(Lists.newArrayList(19),GeneRangeExtractor.extractExonNumber("EGFR exon 19 insertions"));
        assertEquals(Lists.newArrayList(20),GeneRangeExtractor.extractExonNumber("ERBB2 proximal exon 20"));
        assertEquals(Lists.newArrayList(9,11,13,14,17),GeneRangeExtractor.extractExonNumber("KIT mutation in exon 9,11,13,14 or 17"));
        assertEquals(Lists.newArrayList(16,17,18,19),GeneRangeExtractor.extractExonNumber("MET mutation in exon 16-19"));
        assertEquals(Lists.newArrayList(2,3),GeneRangeExtractor.extractExonNumber("Null (Partial deletion of Exons 2 & 3)"));
        assertEquals(Lists.newArrayList(12),GeneRangeExtractor.extractExonNumber("Exon 12 splice site insertion"));

    }

}