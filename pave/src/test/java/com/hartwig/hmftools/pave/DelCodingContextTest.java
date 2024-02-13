package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_DELETION;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createMockGenome;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createNegTranscript;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createPosTranscript;
import static com.hartwig.hmftools.pave.ImpactTestUtils.generateTestBases;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.pave.impact.ImpactClassifier;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

import org.junit.Test;

public class DelCodingContextTest
{
    @Test
    public void testInframeDels()
    {
        MockRefGenome refGenome = createMockGenome();
        String refBases = refGenome.RefGenomeMap.get(CHR_1);

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        TranscriptData transDataPos = createPosTranscript();

        // DELs crossing codons
        int pos = 30;
        String refCodonBases = refBases.substring(pos, pos + 9);

        String ref = refCodonBases.substring(0, 4);
        String alt = refCodonBases.substring(0, 1);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPos);

        assertFalse(impact.codingContext().IsFrameShift);
        assertEquals(refCodonBases.substring(0, 6), impact.proteinContext().RefCodonBases);

        String altCodonBases = alt + refCodonBases.substring(4, 6);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);
        assertEquals(INFRAME_DELETION, impact.topEffect());

        ref = refCodonBases.substring(1, 5);
        alt = refCodonBases.substring(1, 2);
        var = new VariantData(CHR_1, 31, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        assertFalse(impact.codingContext().IsFrameShift);
        assertEquals(refCodonBases.substring(0, 6), impact.proteinContext().RefCodonBases);

        altCodonBases = refCodonBases.substring(0, 1) + alt + refCodonBases.substring(5, 6);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        ref = refCodonBases.substring(2, 6);
        alt = refCodonBases.substring(2, 3);
        var = new VariantData(CHR_1, 32, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        assertFalse(impact.codingContext().IsFrameShift);
        assertEquals(refCodonBases, impact.proteinContext().RefCodonBases);

        altCodonBases = refCodonBases.substring(0, 2) + alt + refCodonBases.substring(6, 9);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);
        assertEquals(INFRAME_DELETION, impact.topEffect());

        // repeat for negative strand
        TranscriptData transDataNeg = createNegTranscript();

        // for codons 80-78 and 77-75
        pos = 72;
        refCodonBases = refBases.substring(pos, pos + 9);

        ref = refCodonBases.substring(0, 4);
        alt = refCodonBases.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        assertFalse(impact.codingContext().IsFrameShift);
        assertEquals(refCodonBases.substring(0, 6), impact.proteinContext().RefCodonBases);

        altCodonBases = alt + refCodonBases.substring(4, 6);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);
        assertEquals(INFRAME_DELETION, impact.topEffect());

        // shift forward one base each time
        ref = refCodonBases.substring(1, 5);
        alt = refCodonBases.substring(1, 2);
        var = new VariantData(CHR_1, ++pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        assertFalse(impact.codingContext().IsFrameShift);
        assertEquals(refCodonBases.substring(0, 6), impact.proteinContext().RefCodonBases);

        altCodonBases = refCodonBases.substring(0, 1) + alt + refCodonBases.substring(5, 6);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        ref = refCodonBases.substring(2, 6);
        alt = refCodonBases.substring(2, 3);
        var = new VariantData(CHR_1, ++pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        assertFalse(impact.codingContext().IsFrameShift);
        assertEquals(refCodonBases, impact.proteinContext().RefCodonBases);

        altCodonBases = refCodonBases.substring(0, 2) + alt + refCodonBases.substring(6, 9);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);
        assertEquals(INFRAME_DELETION, impact.topEffect());
    }

    @Test
    public void testDelsAcrossSplice()
    {
        // DELs crossing exon boundaries
        MockRefGenome refGenome = createMockGenome();
        String refBases = refGenome.RefGenomeMap.get(CHR_1);

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        TranscriptData transDataPos = createPosTranscript();

        int pos = 20;
        String ref = refBases.substring(pos, pos + 4);
        String alt = refBases.substring(pos, pos + 1);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPos);

        assertEquals(CODING, impact.codingContext().CodingType);
        assertEquals(INTRONIC, impact.codingContext().RegionType);
        assertEquals(6, impact.codingContext().CodingBase);
        assertFalse(impact.codingContext().IsFrameShift);

        pos = 19;
        ref = refBases.substring(pos, pos + 4);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        assertTrue(impact.codingContext().IsFrameShift);

        // frameshift since 18 remains while 19-20 are deleted
        pos = 18;
        ref = refBases.substring(pos, pos + 6);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);
        assertTrue(impact.codingContext().IsFrameShift);

        // inframe
        pos = 17;
        ref = refBases.substring(pos, pos + 6);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);
        assertFalse(impact.codingContext().IsFrameShift);

        assertEquals(refBases.substring(15, 21), impact.proteinContext().RefCodonBases);

        String altCodonBases = refBases.substring(15, 18);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);
        assertTrue(impact.effects().contains(INFRAME_DELETION));

        // starting before the next exon
        pos = 28;
        ref = refBases.substring(pos, pos + 7);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);
        assertTrue(impact.codingContext().IsFrameShift);

        pos = 29; // deletes first 3 bases of exon
        ref = refBases.substring(pos, pos + 4);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);
        assertFalse(impact.codingContext().IsFrameShift);
        assertEquals(7, impact.codingContext().CodingBase);
        assertEquals(30, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(33, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals("c.7_9delTCG", impact.codingContext().Hgvs);

        // inframe
        pos = 28;
        ref = refBases.substring(pos, pos + 5);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);
        assertFalse(impact.codingContext().IsFrameShift);

        assertEquals(refBases.substring(30, 36), impact.proteinContext().RefCodonBases);

        altCodonBases = refBases.substring(33, 36);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);
        assertTrue(impact.effects().contains(INFRAME_DELETION));

        // negative strand
        TranscriptData transDataNeg = createNegTranscript();

        pos = 88;
        ref = refBases.substring(pos, pos + 4);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);
        assertTrue(impact.codingContext().IsFrameShift);

        pos = 88;
        ref = refBases.substring(pos, pos + 5);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);
        assertFalse(impact.codingContext().IsFrameShift);

        pos = 78; // deletes last 2 bases of exon, which are first 2 after splice acceptor, index 7 & 8
        ref = refBases.substring(pos, pos + 7);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);
        assertTrue(impact.codingContext().IsFrameShift);
        assertEquals(7, impact.codingContext().CodingBase);
        assertEquals(78, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(80, impact.codingContext().CodingPositionRange[SE_END]);

        pos = 76; // deletes last 3 bases of exon but none of the intron, doesn't span splice
        ref = refBases.substring(pos, pos + 5);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);
        assertTrue(impact.codingContext().IsFrameShift);
        assertEquals(7, impact.codingContext().CodingBase);
        assertEquals(76, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(80, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals("c.7_10delCGAT", impact.codingContext().Hgvs);

        // deletes splice acceptor covering coding and intronic bases, coding bases 7-9
        pos = 77;
        ref = refBases.substring(pos, pos + 5);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);
        assertFalse(impact.codingContext().IsFrameShift);
        assertTrue(impact.codingContext().SpansSpliceJunction);
        assertEquals(3, impact.codingContext().DeletedCodingBases);
        assertEquals(-2, impact.codingContext().NearestExonDistance);
        assertEquals("c.7-1_9delTCGA", impact.codingContext().Hgvs);

        pos = 67;
        ref = refBases.substring(pos, pos + 6);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);
        assertFalse(impact.codingContext().IsFrameShift);
        assertTrue(impact.codingContext().SpansSpliceJunction);
        assertEquals(3, impact.codingContext().DeletedCodingBases);
        assertEquals(3, impact.codingContext().NearestExonDistance);
        assertEquals("c.15_17+2delCGATC", impact.codingContext().Hgvs);

        // stop-codon deletion example
        pos = 12;
        ref = refBases.substring(pos, pos + 6);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);
        // assertTrue(impact.codingContext().IsFrameShift);
        assertEquals(15, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(18, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(3, impact.codingContext().DeletedCodingBases);
        assertEquals(3, impact.codingContext().NearestExonDistance);
        assertTrue(impact.codingContext().SpansCodingEnd);
        assertEquals("c.43_*2delTCGAT", impact.codingContext().Hgvs);
    }

    @Test
    public void testUtrVariants()
    {
        // DELs crossing exon boundaries
        MockRefGenome refGenome = createMockGenome(150);
        String refBases = refGenome.RefGenomeMap.get(CHR_1);

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        TranscriptData transDataNeg = createNegTranscript();

        // DEL in exon before coding begins - 110 = coding-base 6, pos 112 at cb 8, so actually deletes 9-11
        int pos = 112;
        String ref = refBases.substring(pos, pos + 4);
        String alt = ref.substring(0, 1);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataNeg);
        assertEquals(12, impact.codingContext().CodingBase);
        assertEquals(112, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(116, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(3, impact.codingContext().DeletedCodingBases);
        assertEquals("c.-11_-9delGAT", impact.codingContext().Hgvs);

        // 5'UTR delete of exonic bases spanning into next intron - coding base 6 at 110, and 9 at 113, deleted bases are 110-112
        pos = 109;
        ref = refBases.substring(pos, pos + 4);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);
        assertEquals(9, impact.codingContext().CodingBase);
        assertEquals(3, impact.codingContext().DeletedCodingBases);
        assertEquals(110, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(113, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals("c.-8_-6delCGA", impact.codingContext().Hgvs);

        // starting in the next intron
        pos = 107;
        ref = refBases.substring(pos, pos + 6);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);
        assertEquals(9, impact.codingContext().CodingBase);
        assertEquals(3, impact.codingContext().DeletedCodingBases);
        assertEquals(110, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(113, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals("c.-8_-6+2delCGATC", impact.codingContext().Hgvs);

        // spanning into the next intro
        pos = 118;
        ref = refBases.substring(pos, pos + 6);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);
        assertEquals(16, impact.codingContext().CodingBase);
        assertEquals(2, impact.codingContext().DeletedCodingBases);
        assertEquals(118, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(120, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals("c.-16-3_-15delGATCG", impact.codingContext().Hgvs);

        // delete exonic bases just before coding
        TranscriptData transDataPos = createPosTranscript();

        pos = 12;
        ref = refBases.substring(pos, pos + 5);
        alt = refBases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);
        assertEquals(1, impact.codingContext().CodingBase);
        assertTrue(impact.codingContext().SpansCodingStart);
        assertEquals(2, impact.codingContext().DeletedCodingBases);
        assertEquals(15, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(17, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals("c.-2_2delATCG", impact.codingContext().Hgvs);
    }

    @Test
    public void testInframeIndelImpacts()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateTestBases(300);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        int[] exonStarts = { 0, 100, 200 };
        Integer codingStart = Integer.valueOf(10);
        Integer codingEnd = Integer.valueOf(250);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 80, codingStart, codingEnd, false, "");

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        // inframe conservative deletion of 2 codons (ie 3 -> 1)
        int pos = 48; // 3rd base of codon
        String ref = chr1Bases.substring(pos, pos + 7);
        String alt = chr1Bases.substring(pos, pos + 1);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(INFRAME_DELETION, impact.topEffect());
        assertEquals("c.40_45delATCGAT", impact.codingContext().Hgvs);

        // inframe disruptive (spans first and last codon boundaries) deletion of 3 codons (ie 4 -> 1)
        pos = 47; // first base of codon
        ref = chr1Bases.substring(pos, pos + 10);
        alt = chr1Bases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(INFRAME_DELETION, impact.topEffect());
    }

}
