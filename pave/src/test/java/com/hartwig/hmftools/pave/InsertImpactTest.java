package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_INSERTION;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createMockGenome;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createPosTranscript;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class InsertImpactTest
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
        assertEquals(refBases.substring(18, 21), impact.proteinContext().RefCodonBases);

        String altCodonBases = refBases.substring(18, 19) + alt + refBases.substring(20, 21);
        assertEquals(altCodonBases, impact.proteinContext().AltCodonBases);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

    }
}
