package com.hartwig.hmftools.sage.impact;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.test.MockRefGenome.getNextBase;
import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INTRON_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.MISSENSE_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.UPSTREAM_GENE_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.UTR_VARIANT;

import static junit.framework.TestCase.assertEquals;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class VariantImpactTest
{

    @Test
    public void testNonCodingImpacts()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String refBases = generateRandomBases(500);

        refGenome.RefGenomeMap.put(CHR_1, refBases);

        int[] exonStarts = {100, 200, 300, 400};
        Integer codingStart = new Integer(125);
        Integer codingEnd = new Integer(425);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 50, codingStart, codingEnd, false, "");

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        // pre-gene
        int pos = 50;
        VariantData var = createSnv(pos, refBases);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(UPSTREAM_GENE_VARIANT, impact.Consequence);

        // 5' UTR
        pos = 120;
        var = createSnv(pos, refBases);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(UTR_VARIANT, impact.Consequence);

        // 3' UTR
        pos = 440;
        var = createSnv(pos, refBases);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(UTR_VARIANT, impact.Consequence);

        // intronic
        pos = 175;
        var = createSnv(pos, refBases);

        impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(INTRON_VARIANT, impact.Consequence);
    }

    private VariantData createSnv(int position, final String refBases)
    {
        String ref = refBases.substring(position, position + 1);
        String alt = getNextBase(ref);
        return new VariantData(CHR_1, position, ref, alt);
    }

    /*
    @Test
    public void testMissenseImpacts()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateRandomBases(120);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        int[] exonStarts = {0, 20, 40, 60, 80, 100};
        Integer codingStart = new Integer(25);
        Integer codingEnd = new Integer(85);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        int pos = 28;
        String ref = chr1Bases.substring(pos, pos + 1);
        String alt = "G";
        VariantData var = new VariantData(CHR_1, 28, ref, alt);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPosStrand);
        assertEquals(MISSENSE_VARIANT, impact.Consequence);

    }

    private String findMissenseBase()
    {
        return "";
    }
    */
}
