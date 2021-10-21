package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseStrandBases;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SYNONYMOUS;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createMockGenome;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createNegTranscript;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createPosTranscript;
import static com.hartwig.hmftools.pave.ImpactTestUtils.generateAlt;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class MnvCodingContextTest
{
    @Test
    public void testNonExonicPositions()
    {
        MockRefGenome refGenome = createMockGenome();
        String refBases = refGenome.RefGenomeMap.get(CHR_1);

        int[] exonStarts = { 0, 30, 60, 90 };

        TranscriptData transDataPos = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, 35, 65, false, "");

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        // test 5'UTR
        int pos = 22;
        String ref = refBases.substring(pos, pos + 4);
        String alt = generateAlt(ref);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPos);

        // first check general coding context fields
        assertEquals(5, impact.codingContext().CodingBase); // since closer
        assertEquals(UTR_5P, impact.codingContext().CodingType);
        assertEquals(INTRONIC, impact.codingContext().RegionType);
        assertEquals(2, impact.codingContext().ExonRank);
        assertEquals(-8, impact.codingContext().NearestExonDistance);

        pos = 15;
        ref = refBases.substring(pos, pos + 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        // first check general coding context fields
        assertEquals(6, impact.codingContext().CodingBase); // since closer
        assertEquals(1, impact.codingContext().ExonRank);
        assertEquals(5, impact.codingContext().NearestExonDistance);

        // coding
        pos = 45;
        ref = refBases.substring(pos, pos + 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        assertEquals(6, impact.codingContext().CodingBase); // since closer to previous
        assertEquals(CODING, impact.codingContext().CodingType);
        assertEquals(INTRONIC, impact.codingContext().RegionType);
        assertEquals(2, impact.codingContext().ExonRank);
        assertEquals(5, impact.codingContext().NearestExonDistance);

        pos = 51;
        ref = refBases.substring(pos, pos + 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        assertEquals(7, impact.codingContext().CodingBase); // since closer to previous
        assertEquals(CODING, impact.codingContext().CodingType);
        assertEquals(INTRONIC, impact.codingContext().RegionType);
        assertEquals(3, impact.codingContext().ExonRank);
        assertEquals(-9, impact.codingContext().NearestExonDistance);

        // 3'UTR
        pos = 72;
        ref = refBases.substring(pos, pos + 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        assertEquals(5, impact.codingContext().CodingBase); // since closer to previous
        assertEquals(UTR_3P, impact.codingContext().CodingType);
        assertEquals(INTRONIC, impact.codingContext().RegionType);
        assertEquals(3, impact.codingContext().ExonRank);
        assertEquals(2, impact.codingContext().NearestExonDistance);

        pos = 81;
        ref = refBases.substring(pos, pos + 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        assertEquals(6, impact.codingContext().CodingBase); // since closer to previous
        assertEquals(4, impact.codingContext().ExonRank);
        assertEquals(-9, impact.codingContext().NearestExonDistance);

        // repeat for negative strand
        TranscriptData transDataNeg = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, exonStarts, 10, 35, 65, false, "");

        // test 5'UTR
        pos = 85;
        ref = refBases.substring(pos, pos + 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        // first check general coding context fields
        assertEquals(6, impact.codingContext().CodingBase); // since closert to upstream exon
        assertEquals(UTR_5P, impact.codingContext().CodingType);
        assertEquals(INTRONIC, impact.codingContext().RegionType);
        assertEquals(1, impact.codingContext().ExonRank);
        assertEquals(2, impact.codingContext().NearestExonDistance);

        pos = 73;
        ref = refBases.substring(pos, pos + 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        // first check general coding context fields
        assertEquals(5, impact.codingContext().CodingBase);
        assertEquals(2, impact.codingContext().ExonRank);
        assertEquals(-6, impact.codingContext().NearestExonDistance);

        // coding
        pos = 51;
        ref = refBases.substring(pos, pos + 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        assertEquals(6, impact.codingContext().CodingBase);
        assertEquals(CODING, impact.codingContext().CodingType);
        assertEquals(INTRONIC, impact.codingContext().RegionType);
        assertEquals(2, impact.codingContext().ExonRank);
        assertEquals(6, impact.codingContext().NearestExonDistance);

        pos = 45;
        ref = refBases.substring(pos, pos + 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        assertEquals(7, impact.codingContext().CodingBase); // since closer to previous
        assertEquals(CODING, impact.codingContext().CodingType);
        assertEquals(INTRONIC, impact.codingContext().RegionType);
        assertEquals(3, impact.codingContext().ExonRank);
        assertEquals(-8, impact.codingContext().NearestExonDistance);

        // 3'UTR
        pos = 22;
        ref = refBases.substring(pos, pos + 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        assertEquals(5, impact.codingContext().CodingBase); // since closer to previous
        assertEquals(UTR_3P, impact.codingContext().CodingType);
        assertEquals(INTRONIC, impact.codingContext().RegionType);
        assertEquals(3, impact.codingContext().ExonRank);
        assertEquals(5, impact.codingContext().NearestExonDistance);

        pos = 15;
        ref = refBases.substring(pos, pos + 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        assertEquals(6, impact.codingContext().CodingBase); // since closer to previous
        assertEquals(4, impact.codingContext().ExonRank);
        assertEquals(-8, impact.codingContext().NearestExonDistance);
    }

    @Test
    public void testMnvBasic()
    {
        MockRefGenome refGenome = createMockGenome();
        String refBases = refGenome.RefGenomeMap.get(CHR_1);

        TranscriptData transDataPos = createPosTranscript();

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        // MNV spanning 2 codons at all positions - codons 30-32 and 33-35
        int pos = 30;
        String refCodonBases = refBases.substring(pos, pos + 6);

        String ref = refCodonBases.substring(0, 4);
        String alt = generateAlt(ref);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPos);

        // first check general coding context fields
        assertEquals(7, impact.codingContext().CodingBase);
        assertEquals(30, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(33, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(PHASE_1, impact.codingContext().UpstreamPhase);
        assertEquals(2, impact.codingContext().ExonRank);
        assertEquals(CODING, impact.codingContext().CodingType);
        assertEquals(EXONIC, impact.codingContext().RegionType);

        // then the protein context
        assertTrue(impact.proteinContext() != null);
        assertEquals(3, impact.proteinContext().CodonIndex);

        String altCodonBases = alt + refCodonBases.substring(4, 6);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        // shift up one base
        pos = 31;
        ref = refCodonBases.substring(1, 5);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);
        impact = classifier.classifyVariant(var, transDataPos);

        assertEquals(8, impact.codingContext().CodingBase);
        assertEquals(31, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(34, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(PHASE_2, impact.codingContext().UpstreamPhase);
        assertEquals(3, impact.proteinContext().CodonIndex);

        altCodonBases = refCodonBases.substring(0, 1) + alt + refCodonBases.substring(5, 6);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        pos = 32;
        ref = refCodonBases.substring(2, 6);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);
        impact = classifier.classifyVariant(var, transDataPos);

        assertEquals(9, impact.codingContext().CodingBase);
        assertEquals(32, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(35, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(PHASE_0, impact.codingContext().UpstreamPhase);
        assertEquals(3, impact.proteinContext().CodonIndex);

        altCodonBases = refCodonBases.substring(0, 2) + alt;
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        // again but on negative strand
        TranscriptData transDataNeg = createNegTranscript();

        // MNV spanning 2 codons at all positions - codons 80-78 and 77-75
        pos = 75;
        refCodonBases = refBases.substring(pos, pos + 6);

        ref = refCodonBases.substring(2, 6);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, 77, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        // first check general coding context fields
        assertEquals(7, impact.codingContext().CodingBase);
        assertEquals(77, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(80, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(PHASE_1, impact.codingContext().UpstreamPhase);
        assertEquals(2, impact.codingContext().ExonRank);
        assertEquals(CODING, impact.codingContext().CodingType);
        assertEquals(EXONIC, impact.codingContext().RegionType);

        // then the protein context
        assertTrue(impact.proteinContext() != null);
        assertEquals(3, impact.proteinContext().CodonIndex);

        altCodonBases = refCodonBases.substring(0, 2) + alt;
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        // shift along 1 each time
        ref = refCodonBases.substring(1, 5);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, 76, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        // first check general coding context fields
        assertEquals(8, impact.codingContext().CodingBase);
        assertEquals(76, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(79, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(PHASE_2, impact.codingContext().UpstreamPhase);
        assertEquals(3, impact.proteinContext().CodonIndex);

        altCodonBases = refCodonBases.substring(0, 1) + alt + refCodonBases.substring(5, 6);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        // and again
        ref = refCodonBases.substring(0, 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, 75, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        // first check general coding context fields
        assertEquals(9, impact.codingContext().CodingBase);
        assertEquals(75, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(78, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(PHASE_0, impact.codingContext().UpstreamPhase);
        assertEquals(3, impact.proteinContext().CodonIndex);

        altCodonBases = alt + refCodonBases.substring(4, 6);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);
    }

    @Test
    public void testMnvAcrossSplice()
    {
        MockRefGenome refGenome = createMockGenome();

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        // MNV spanning exon boundary - shortens the effect
        String refBases = refGenome.RefGenomeMap.get(CHR_1);

        TranscriptData transDataPos = createPosTranscript();

        // firstly a 2-lots across the exon
        String refCodonBases = refBases.substring(15, 21);

        int pos = 20;
        String ref = refBases.substring(pos, pos + 2);
        String alt = generateAlt(ref);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPos);

        // first check general coding context fields
        assertEquals(6, impact.codingContext().CodingBase);
        assertEquals(20, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(20, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(PHASE_0, impact.codingContext().UpstreamPhase);
        assertEquals(1, impact.codingContext().ExonRank);
        assertEquals(CODING, impact.codingContext().CodingType);
        assertTrue(impact.codingContext().SpansSpliceJunction);
        assertEquals(EXONIC, impact.codingContext().RegionType);

        assertTrue(impact.proteinContext() != null);
        assertEquals(2, impact.proteinContext().CodonIndex);
        assertEquals(refCodonBases.substring(3, 6), impact.proteinContext().RefCodonBases);

        String altCodonBases = refCodonBases.substring(3, 5) + alt.substring(0, 1);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        // same again but spanning 2 exons with the codons: 36-38, 39-50, 51-53
        refCodonBases = refBases.substring(36, 41) + refBases.substring(50, 54);

        pos = 38;
        ref = refBases.substring(pos, pos + 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        // first check general coding context fields
        assertEquals(15, impact.codingContext().CodingBase);
        assertEquals(38, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(40, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(PHASE_0, impact.codingContext().UpstreamPhase);
        assertEquals(2, impact.codingContext().ExonRank);
        assertEquals(CODING, impact.codingContext().CodingType);
        assertEquals(EXONIC, impact.codingContext().RegionType);
        assertEquals("c.15_17+1delTCGAinsCGAT", impact.codingContext().Hgvs);
        assertTrue(impact.proteinContext() != null);
        assertEquals(5, impact.proteinContext().CodonIndex);
        assertEquals(refBases.substring(36, 41) + refBases.substring(50, 51), impact.proteinContext().RefCodonBases);

        altCodonBases = refBases.substring(36, 38) + alt.substring(0, 3) + refBases.substring(50, 51);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        // and starting before the next exon
        pos = 48;
        ref = refBases.substring(pos, pos + 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        // first check general coding context fields
        assertEquals(18, impact.codingContext().CodingBase);
        assertEquals(50, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(51, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(PHASE_0, impact.codingContext().UpstreamPhase);
        assertEquals(3, impact.codingContext().ExonRank);
        assertEquals(CODING, impact.codingContext().CodingType);
        assertEquals(EXONIC, impact.codingContext().RegionType);
        assertEquals("c.18-2_19delGATCinsATCG", impact.codingContext().Hgvs);

        assertTrue(impact.proteinContext() != null);
        assertEquals(6, impact.proteinContext().CodonIndex);
        assertEquals(refBases.substring(39, 41) + refBases.substring(50, 54), impact.proteinContext().RefCodonBases);

        altCodonBases = refBases.substring(39, 41) + alt.substring(2, 4) + refBases.substring(52, 54);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        // for negative strand, is the phase set correctly even though the variant starts in the upstream intron?

        TranscriptData transDataNeg = createNegTranscript();

        // and starting at the end of the last exon
        pos = 88;
        ref = refBases.substring(pos, pos + 3);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        // first check general coding context fields
        assertEquals(6, impact.codingContext().CodingBase);
        assertEquals(90, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(90, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(PHASE_0, impact.codingContext().UpstreamPhase);
        assertEquals(1, impact.codingContext().ExonRank);
        assertTrue(impact.codingContext().SpansSpliceJunction);
        assertEquals(CODING, impact.codingContext().CodingType);
        assertEquals(EXONIC, impact.codingContext().RegionType);

        assertTrue(impact.proteinContext() != null);
        assertEquals(2, impact.proteinContext().CodonIndex);
        assertEquals(refBases.substring(90, 93), impact.proteinContext().RefCodonBases);

        altCodonBases = alt.substring(2, 3) + refBases.substring(91, 93);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        // spanning exons 2 and 3 codons: 74-72, 71-60, 59-57
        pos = 68;
        ref = refBases.substring(pos, pos + 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        // first check general coding context fields
        assertEquals(16, impact.codingContext().CodingBase);
        assertEquals(70, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(71, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(PHASE_1, impact.codingContext().UpstreamPhase);
        assertEquals(2, impact.codingContext().ExonRank);
        assertTrue(impact.codingContext().SpansSpliceJunction);

        assertTrue(impact.proteinContext() != null);
        assertEquals(6, impact.proteinContext().CodonIndex);
        assertEquals(refBases.substring(60, 61) + refBases.substring(70, 72), impact.proteinContext().RefCodonBases);

        altCodonBases = refBases.substring(60, 61) + alt.substring(2, 4);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);

        pos = 59;
        ref = refBases.substring(pos, pos + 4);
        alt = generateAlt(ref);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        // first check general coding context fields
        assertEquals(18, impact.codingContext().CodingBase);
        assertEquals(59, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(60, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(PHASE_0, impact.codingContext().UpstreamPhase);
        assertEquals(3, impact.codingContext().ExonRank);
        assertTrue(impact.codingContext().SpansSpliceJunction);

        assertTrue(impact.proteinContext() != null);
        assertEquals(6, impact.proteinContext().CodonIndex);
        assertEquals(refBases.substring(57, 61) + refBases.substring(70, 72), impact.proteinContext().RefCodonBases);

        altCodonBases = refBases.substring(57, 59) + alt.substring(0, 2) + refBases.substring(70, 72);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);
    }

}
