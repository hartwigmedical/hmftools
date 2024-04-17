package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_INSERTION;
import static com.hartwig.hmftools.pave.impact.HgvsProtein.HGVS_SPLICE_UNKNOWN;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createMockGenome;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createNegTranscript;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createPosTranscript;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.pave.impact.ImpactClassifier;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

import org.junit.Test;

public class InsertCodingContextTest
{
    @Test
    public void testInserts()
    {
        MockRefGenome refGenome = createMockGenome();
        String refBases = refGenome.RefGenomeMap.get(CHR_1);

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        TranscriptData transDataPos = createPosTranscript();

        // insert at last base of exon does nothing
        int pos = 20;
        String ref = refBases.substring(20, 21);
        String alt = ref + "AAA";
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPos);

        assertEquals(CODING, impact.codingContext().CodingType);
        assertEquals(INTRONIC, impact.codingContext().RegionType);
        assertEquals(6, impact.codingContext().CodingBase);
        assertEquals(HGVS_SPLICE_UNKNOWN, impact.proteinContext().Hgvs);

        // frameshift
        pos = 19;
        ref = refBases.substring(19, 20);
        alt = ref + "AAAA";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        assertEquals(5, impact.codingContext().CodingBase);
        assertEquals(19, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(19, impact.codingContext().CodingPositionRange[SE_END]);
        assertTrue(impact.codingContext().IsFrameShift);

        // inframe insert
        pos = 19;
        ref = refBases.substring(19, 20);
        alt = ref + "AAA";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        assertEquals(5, impact.codingContext().CodingBase);
        assertEquals(19, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(19, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(PHASE_2, impact.codingContext().UpstreamPhase);
        assertEquals(1, impact.codingContext().ExonRank);
        assertEquals(CODING, impact.codingContext().CodingType);
        assertEquals(EXONIC, impact.codingContext().RegionType);

        assertFalse(impact.codingContext().IsFrameShift);
        assertEquals(refBases.substring(15, 21) + refBases.substring(30, 33), impact.proteinContext().RefCodonBases);

        String altCodonBases = refBases.substring(15, 19) + alt + refBases.substring(20, 21) + refBases.substring(30, 33);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        // inframe insert at exon start
        pos = 50;
        ref = refBases.substring(50, 51);
        alt = ref + "AAA";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        assertEquals(18, impact.codingContext().CodingBase);
        assertEquals(50, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(50, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(3, impact.codingContext().ExonRank);

        assertFalse(impact.codingContext().IsFrameShift);
        assertEquals(refBases.substring(39, 41) + refBases.substring(50, 54), impact.proteinContext().RefCodonBases);

        altCodonBases = refBases.substring(39, 41) + alt + refBases.substring(51, 54);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        // insert just before exon start
        pos = 29;
        ref = refBases.substring(29, 30);
        alt = ref + "A";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);

        assertEquals(7, impact.codingContext().CodingBase);
        assertEquals(30, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(30, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(2, impact.codingContext().ExonRank);

        assertTrue(impact.codingContext().IsFrameShift);
        assertTrue(impact.inSpliceRegion());
        assertEquals(FRAMESHIFT, impact.topEffect());

        // repeat for negative strand
        TranscriptData transDataNeg = createNegTranscript();

        // at exon boundary - inserted bases go into coding region
        pos = 80;
        ref = refBases.substring(80, 81);
        alt = ref + "CGG";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        assertEquals(CODING, impact.codingContext().CodingType);
        assertEquals(EXONIC, impact.codingContext().RegionType);
        assertEquals(7, impact.codingContext().CodingBase);
        assertEquals(80, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(80, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(4, impact.codingContext().ExonRank);

        assertTrue(impact.inSpliceRegion());
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        // frameshift
        pos = 91;
        ref = refBases.substring(91, 92);
        alt = ref + "AAAA";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        assertEquals(4, impact.codingContext().CodingBase);
        assertEquals(92, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(92, impact.codingContext().CodingPositionRange[SE_END]);
        assertTrue(impact.codingContext().IsFrameShift);

        // inframe insert
        pos = 91;
        ref = refBases.substring(91, 92);
        alt = ref + "GGG";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        assertEquals(4, impact.codingContext().CodingBase);
        assertEquals(92, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(92, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(PHASE_1, impact.codingContext().UpstreamPhase);
        assertEquals(3, impact.codingContext().ExonRank);
        assertEquals(CODING, impact.codingContext().CodingType);
        assertEquals(EXONIC, impact.codingContext().RegionType);

        assertFalse(impact.codingContext().IsFrameShift);
        assertEquals(refBases.substring(78, 81) + refBases.substring(90, 96), impact.proteinContext().RefCodonBases);

        altCodonBases = refBases.substring(78, 81) + refBases.substring(90, 91) + alt + refBases.substring(92, 96);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        // inframe insert at exon end
        pos = 70;
        ref = refBases.substring(70, 71);
        alt = ref + "AAA";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        assertEquals(16, impact.codingContext().CodingBase);
        assertEquals(71, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(71, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(4, impact.codingContext().ExonRank);

        assertFalse(impact.codingContext().IsFrameShift);
        assertEquals(refBases.substring(57, 61) + refBases.substring(70, 75), impact.proteinContext().RefCodonBases);

        altCodonBases = refBases.substring(57, 61) + alt + refBases.substring(71, 75);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        // test at phase 0
        pos = 92;
        ref = refBases.substring(92, 93);
        alt = ref + "GGG";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);

        assertEquals(3, impact.codingContext().CodingBase);
        assertEquals(93, impact.codingContext().CodingPositionRange[SE_START]);
        assertEquals(93, impact.codingContext().CodingPositionRange[SE_END]);
        assertEquals(PHASE_0, impact.codingContext().UpstreamPhase);
        assertEquals(3, impact.codingContext().ExonRank);

        assertFalse(impact.codingContext().IsFrameShift);
        assertEquals(refBases.substring(90, 96), impact.proteinContext().RefCodonBases);

        altCodonBases = refBases.substring(90, 93) + alt.substring(1) + refBases.substring(93, 96);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);
        assertEquals(INFRAME_INSERTION, impact.topEffect());
    }
}
