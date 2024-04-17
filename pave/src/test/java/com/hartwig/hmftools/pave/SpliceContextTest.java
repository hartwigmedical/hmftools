package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.pave.ImpactTestUtils.createMockGenome;
import static com.hartwig.hmftools.pave.ImpactTestUtils.generateAlt;
import static com.hartwig.hmftools.pave.impact.SpliceClassifier.determineIndelSpliceImpact;
import static com.hartwig.hmftools.pave.impact.SpliceImpactType.BASE_SHIFT;
import static com.hartwig.hmftools.pave.impact.SpliceImpactType.OUTSIDE_RANGE;
import static com.hartwig.hmftools.pave.impact.SpliceImpactType.REGION_DELETED;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.pave.impact.ImpactClassifier;
import com.hartwig.hmftools.pave.impact.SpliceImpactType;

import org.junit.Test;

public class SpliceContextTest
{
    private MockRefGenome mRefGenome = createMockGenome(41);
    private String mRefBases = mRefGenome.RefGenomeMap.get(CHR_1);

    private boolean refAltSpliceBasesMatch(
            final VariantData variant, final RefGenomeInterface refGenome, final ExonData exon, boolean posStrand)
    {
        SpliceImpactType type = determineIndelSpliceImpact(variant, refGenome, exon, posStrand);
        return type != BASE_SHIFT;
    }

    @Test
    public void testDeleteSpliceBasesPosStrand()
    {
        TranscriptData transDataPos = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] { 20 }, 10, 20, 30, false, "");

        ExonData exon = transDataPos.exons().get(0);

        // test deleting the whole splice region
        int pos = 15;
        String ref = mRefBases.substring(pos, pos + 6);
        String alt = ref.substring(0, 1);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        assertEquals(REGION_DELETED, determineIndelSpliceImpact(var, mRefGenome, exon, transDataPos.posStrand()));

        // and on the donor side
        pos = 29;
        ref = mRefBases.substring(pos, pos + 7);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertEquals(REGION_DELETED, determineIndelSpliceImpact(var, mRefGenome, exon, transDataPos.posStrand()));

        // deletes outside region are non-disruptive
        pos = 15;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        pos = 19;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        pos = 28;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        pos = 35;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        // acceptor positions - modify A3 only but doesn't become a G so is ok
        pos = 17;
        ref = mRefBases.substring(pos, pos + 1);
        alt = "A";
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        // acceptor positions - delete A3 and cannot pull exonic bases, so remains disruptive
        pos = 15;
        ref = mRefBases.substring(pos, pos + 3);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        // deletes A3 and A2
        pos = 16;
        ref = mRefBases.substring(pos, pos + 3);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        // now with A3 becoming disruptive
        pos = 16;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        // delete A2 only
        pos = 17;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        // delete A2 and A1
        pos = 17;
        ref = mRefBases.substring(pos, pos + 3);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        // donor bases - starting with delete of D-1 only
        pos = 29;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        pos = 30;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        pos = 34;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));
    }

    @Test
    public void testDeleteSpliceBasesNegStrand()
    {
        TranscriptData transDataNeg = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, new int[] { 20 }, 10, 20, 30, false, "");

        ExonData exon = transDataNeg.exons().get(0);

        // deletes outside region are non-disruptive
        int pos = 12;
        String ref = mRefBases.substring(pos, pos + 3);
        String alt = ref.substring(0, 1);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        pos = 20;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        pos = 29;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        pos = 33;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        // acceptor positions - delete A3 only but doesn't become a G so is ok
        pos = 32;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        // deletes A3 and A2
        pos = 31;
        ref = mRefBases.substring(pos, pos + 3);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        // now with A3 becoming disruptive
        pos = 32;
        ref = mRefBases.substring(pos, pos + 3);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        // delete A2 only
        pos = 31;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        // delete A2 and A1
        pos = 30;
        ref = mRefBases.substring(pos, pos + 3);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        // donor bases - starting with delete of D-1 only
        pos = 19;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        // D1
        pos = 18;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        // D4
        pos = 15;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        // D5
        pos = 14;
        ref = mRefBases.substring(pos, pos + 2);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));
    }

    @Test
    public void testInsertSpliceBasesPosStrand()
    {
        TranscriptData transDataPos = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] { 20 }, 10, 20, 30, false, "");

        ExonData exon = transDataPos.exons().get(0);

        int pos = 17;
        String ref = mRefBases.substring(pos, pos + 1);
        String alt = ref + generateAlt(mRefBases.substring(pos + 1, pos + 2));
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        ref = mRefBases.substring(pos, pos + 1);
        alt = mRefBases.substring(pos, pos + 2);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        // now insert bases that match the A1-3 bases - these will be considered disruptive before the exon
        ref = mRefBases.substring(pos, pos + 1);
        alt = mRefBases.substring(pos, pos + 3);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        // starting at A2, inserting matching bases
        pos = 18;
        ref = mRefBases.substring(pos, pos + 1);
        alt = mRefBases.substring(17, 20);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        // starting at A1 - makes no difference since inserts into exon
        pos = 19;
        ref = mRefBases.substring(pos, pos + 1);
        alt = "AAAA";
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        // test D-1 through D5 still on positive strand
        pos = 30;
        ref = mRefBases.substring(pos, pos + 1);
        alt = mRefBases.substring(pos, pos + 2);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        pos = 30;
        ref = mRefBases.substring(pos, pos + 1);
        alt = mRefBases.substring(pos, pos + 5);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        pos = 30;
        ref = mRefBases.substring(pos, pos + 1);
        alt = mRefBases.substring(pos, pos + 6);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        // from D2
        pos = 32;
        ref = mRefBases.substring(pos, pos + 1);
        alt = mRefBases.substring(pos, pos + 6);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        pos = 32;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + generateAlt(mRefBases.substring(pos + 1, pos + 4));
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        // from D4
        pos = 34;
        ref = mRefBases.substring(pos, pos + 1);
        alt = mRefBases.substring(pos, pos + 2);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));

        pos = 34;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + generateAlt(mRefBases.substring(pos + 1, pos + 2));
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataPos.posStrand()));
    }

    @Test
    public void testInsertSpliceBasesNegStrand()
    {
        // repeat for negative strand
        TranscriptData transDataNeg = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, new int[] {20}, 10, 20, 30, false, "");

        ExonData exon = transDataNeg.exons().get(0);

        // inserts outside the range are non-disruptive
        int pos = 33;
        String ref = mRefBases.substring(pos, pos + 1);
        String alt = ref + "AAAA";
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        pos = 29;
        ref = mRefBases.substring(pos, pos + 1);
        alt = alt = ref + "AAAA";
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        pos = 14;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + "AAAA";
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        pos = 20;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + "AAAA";
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));


        // first the splice acceptors
        pos = 30;
        ref = mRefBases.substring(pos, pos + 1);
        alt = mRefBases.substring(pos, pos + 4);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertEquals(OUTSIDE_RANGE, determineIndelSpliceImpact(var, mRefGenome, exon, transDataNeg.posStrand()));

        alt = ref + generateAlt(mRefBases.substring(pos + 1, pos + 4));
        var = new VariantData(CHR_1, pos, ref, alt);

        assertEquals(OUTSIDE_RANGE, determineIndelSpliceImpact(var, mRefGenome, exon, transDataNeg.posStrand()));

        // from A2 - first with A3 not changing to the required C
        pos = 32;
        ref = mRefBases.substring(pos, pos + 1);
        alt = mRefBases.substring(pos, pos + 4);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        alt = ref + "CCCC";
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        // then the donors

        // from D5
        pos = 15;
        ref = mRefBases.substring(pos, pos + 1);
        alt = mRefBases.substring(pos, pos + 4);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        // ok if covers the full range with homologous bases
        alt = mRefBases.substring(pos, pos + 6);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertTrue(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        alt = ref + "A";
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));

        // D1
        pos = 19;
        ref = mRefBases.substring(pos, pos + 1);
        alt = mRefBases.substring(pos, pos + 4);
        var = new VariantData(CHR_1, pos, ref, alt);

        assertFalse(refAltSpliceBasesMatch(var, mRefGenome, exon, transDataNeg.posStrand()));
    }
}
