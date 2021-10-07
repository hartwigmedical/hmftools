package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_DELETION;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createMockGenome;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createNegTranscript;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createPosTranscript;
import static com.hartwig.hmftools.pave.ImpactTestUtils.generateAlt;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class DelImpactTest
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

        assertEquals(3, impact.codingContext().DeletedCodingBases);
        assertEquals(refCodonBases.substring(0, 6), impact.proteinContext().RefCodonBases);

        String altCodonBases = alt + refCodonBases.substring(4, 6);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        ref = refCodonBases.substring(1, 5);
        alt = refCodonBases.substring(1, 2);
        var = new VariantData(CHR_1, 31, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        assertEquals(3, impact.codingContext().DeletedCodingBases);
        assertEquals(refCodonBases.substring(0, 6), impact.proteinContext().RefCodonBases);

        altCodonBases = refCodonBases.substring(0, 1) + alt + refCodonBases.substring(5, 6);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        ref = refCodonBases.substring(2, 6);
        alt = refCodonBases.substring(2, 3);
        var = new VariantData(CHR_1, 32, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        assertEquals(3, impact.codingContext().DeletedCodingBases);
        assertEquals(refCodonBases, impact.proteinContext().RefCodonBases);

        altCodonBases = refCodonBases.substring(0, 2) + alt + refCodonBases.substring(6, 9);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        // repeat for negative strand
        TranscriptData transDataNeg = createNegTranscript();

        // for codons 80-78 and 77-75
        pos = 72;
        refCodonBases = refBases.substring(pos, pos + 9);

        ref = refCodonBases.substring(0, 4);
        alt = refCodonBases.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        assertEquals(3, impact.codingContext().DeletedCodingBases);
        assertEquals(refCodonBases.substring(0, 6), impact.proteinContext().RefCodonBases);

        altCodonBases = alt + refCodonBases.substring(4, 6);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        // shift forward one base each time
        ref = refCodonBases.substring(1, 5);
        alt = refCodonBases.substring(1, 2);
        var = new VariantData(CHR_1, ++pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        assertEquals(3, impact.codingContext().DeletedCodingBases);
        assertEquals(refCodonBases.substring(0, 6), impact.proteinContext().RefCodonBases);

        altCodonBases = refCodonBases.substring(0, 1) + alt + refCodonBases.substring(5, 6);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        ref = refCodonBases.substring(2, 6);
        alt = refCodonBases.substring(2, 3);
        var = new VariantData(CHR_1, ++pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        assertEquals(3, impact.codingContext().DeletedCodingBases);
        assertEquals(refCodonBases, impact.proteinContext().RefCodonBases);

        altCodonBases = refCodonBases.substring(0, 2) + alt + refCodonBases.substring(6, 9);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        // exon boundary tests
    }

    @Test
    public void testInframeIndelImpacts()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateRandomBases(300);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        int[] exonStarts = { 0, 100, 200 };
        Integer codingStart = new Integer(10);
        Integer codingEnd = new Integer(250);

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

        // inframe disruptive (spans first and last codon boundaries) deletion of 3 codons (ie 4 -> 1)
        pos = 47; // first base of codon
        ref = chr1Bases.substring(pos, pos + 10);
        alt = chr1Bases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(INFRAME_DELETION, impact.topEffect());

    }

}
