package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_INSERTION;
import static com.hartwig.hmftools.pave.ImpactTestUtils.generateTestBases;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;

import static junit.framework.TestCase.assertEquals;

import com.hartwig.hmftools.common.fusion.FusionCommon;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.pave.impact.CodingContext;
import com.hartwig.hmftools.pave.impact.HgvsCoding;
import com.hartwig.hmftools.pave.impact.ImpactClassifier;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

import org.junit.Test;

public class HgvsCodingTest
{
    @Test
    public void testPointMutations()
    {
        int pos = 100;
        String ref = "A";
        String alt = "C";
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        CodingContext codingContext = new CodingContext();
        codingContext.Strand = POS_STRAND;

        codingContext.RegionType = INTRONIC;
        codingContext.CodingType = UTR_5P;
        codingContext.CodingBase = 10;
        codingContext.NearestExonDistance = -50;

        assertEquals("c.-10-50A>C", HgvsCoding.generate(var, codingContext));

        codingContext.RegionType = INTRONIC;
        codingContext.CodingType = UTR_3P;
        codingContext.CodingBase = 20;
        codingContext.NearestExonDistance = 10;

        assertEquals("c.*20+10A>C", HgvsCoding.generate(var, codingContext));

        codingContext.RegionType = INTRONIC;
        codingContext.CodingType = CODING;
        codingContext.CodingBase = 10;
        codingContext.NearestExonDistance = -50;

        assertEquals("c.10-50A>C", HgvsCoding.generate(var, codingContext));

        codingContext.NearestExonDistance = 50;

        assertEquals("c.10+50A>C", HgvsCoding.generate(var, codingContext));

        codingContext.RegionType = EXONIC;
        assertEquals("c.10A>C", HgvsCoding.generate(var, codingContext));
    }

    @Test
    public void testCodingDuplication()
    {
        // pos codons: M 20-22, A 23-25, C 26-28, D 29-31, L 32-34, L 35-37, G 38-40, H 41-43, E 44-46, stopX 47-49
        // amino acids:           M  A  C  D  L  L  G  H  E  X
        // exon 1:                20        30
        // position:              012345678901234
        String refCodingBases1 = "ATGGCTTGTGACTTA";

        // intron:                   40
        // position             567890123456789
        // bases into intron:   12345   7654321
        String intronicBases = "TTTTGGGGCCCCAAA";

        // exon 2:                50        60
        // position:              012345678901234
        String refCodingBases2 = "TTAGGACACGAGTAA";

        MockRefGenome refGenome = new MockRefGenome();

        String refBases = generateTestBases(20) + refCodingBases1 + intronicBases + refCodingBases2 + generateTestBases(20);
        refGenome.RefGenomeMap.put(CHR_1, refBases);
        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        TranscriptData posTrans = createTransExons(
                GENE_ID_1, TRANS_ID_1, FusionCommon.POS_STRAND, new int[] {10, 50}, 24,
                20, 64, false, "");

        // duplication of a codon
        int pos = 25;
        String ref = refBases.substring(pos, pos + 1);
        String alt = refBases.substring(pos, pos + 4);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        String altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        VariantTransImpact impact = classifier.classifyVariant(var, posTrans);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        assertEquals("c.7_9dupTGT", impact.codingContext().Hgvs);
        assertEquals("p.Cys3dup", impact.proteinContext().Hgvs);

        // duplication of 1 intronic base
        pos = 40;
        ref = refBases.substring(pos, pos + 1);
        alt = refBases.substring(pos, pos + 2);
        var = new VariantData(CHR_1, pos, ref, alt);

        altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        impact = classifier.classifyVariant(var, posTrans);
        assertEquals(VariantEffect.INTRONIC, impact.topEffect());
        assertEquals(6, impact.codingContext().NearestExonDistance);
        assertEquals("c.15+7dupG", impact.codingContext().Hgvs);

        // duplication of 2 intronic bases
        pos = 40;
        ref = refBases.substring(pos, pos + 1);
        alt = refBases.substring(pos, pos + 3);
        var = new VariantData(CHR_1, pos, ref, alt);

        altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        impact = classifier.classifyVariant(var, posTrans);
        assertEquals(VariantEffect.INTRONIC, impact.topEffect());
        assertEquals(6, impact.codingContext().NearestExonDistance);
        assertEquals("c.15+7_15+8dupGG", impact.codingContext().Hgvs);

        // test nearer to the next exon
        pos = 47;
        ref = refBases.substring(pos, pos + 1);
        alt = refBases.substring(pos, pos + 3);
        var = new VariantData(CHR_1, pos, ref, alt);

        altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        impact = classifier.classifyVariant(var, posTrans);
        assertEquals(VariantEffect.SPLICE_ACCEPTOR, impact.topEffect());

        assertEquals(16, impact.codingContext().CodingBase);
        assertEquals(-3, impact.codingContext().NearestExonDistance);

        // the base positions 48-49 are duplicated, so these are intronic bases -1 and -2
        assertEquals("c.16-2_16-1dupAA", impact.codingContext().Hgvs);

        // also for a single duplicated base
        alt = refBases.substring(pos, pos + 2);
        var = new VariantData(CHR_1, pos, ref, alt);

        altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);
        impact = classifier.classifyVariant(var, posTrans);
        assertEquals("c.16-2dupA", impact.codingContext().Hgvs);

        // again on the negative strand
        TranscriptData negTrans = createTransExons(
                GENE_ID_1, TRANS_ID_1, FusionCommon.NEG_STRAND, new int[] {10, 50}, 24,
                20, 64, false, "");

        // coding goes from 64 -> 50 then 34 -> 20

        // first an exonic dup
        pos = 58;
        ref = refBases.substring(pos, pos + 1);
        alt = refBases.substring(pos, pos + 4);
        var = new VariantData(CHR_1, pos, ref, alt);

        altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        impact = classifier.classifyVariant(var, negTrans);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        assertEquals("c.4_6dupCTC", impact.codingContext().Hgvs);

        // dup in intron following an exon
        pos = 43;
        ref = refBases.substring(pos, pos + 1);
        alt = refBases.substring(pos, pos + 4);
        var = new VariantData(CHR_1, pos, ref, alt);

        altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        impact = classifier.classifyVariant(var, negTrans);
        assertEquals(VariantEffect.INTRONIC, impact.topEffect());
        assertEquals(15, impact.codingContext().CodingBase);
        assertEquals(6, impact.codingContext().NearestExonDistance);

        // bases at pos 44-46 are duplicated, these are intronic bases 6-4
        assertEquals("c.15+4_15+6dupGGG", impact.codingContext().Hgvs);

        // single-base dup in prior intron
        pos = 35;
        ref = refBases.substring(pos, pos + 1);
        alt = refBases.substring(pos, pos + 1) + refBases.substring(pos + 1, pos + 2); // duplicating the next base
        var = new VariantData(CHR_1, pos, ref, alt);

        altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        impact = classifier.classifyVariant(var, negTrans);

        assertEquals(16, impact.codingContext().CodingBase);
        assertEquals(-2 , impact.codingContext().NearestExonDistance); // nearest exon ends at 34, using upstream pos of 37

        assertEquals("c.16-2dupA", impact.codingContext().Hgvs);

        pos = 36;
        ref = refBases.substring(pos, pos + 1);
        alt = ref + refBases.substring(pos + 1, pos + 2) + refBases.substring(pos + 1, pos + 2);
        var = new VariantData(CHR_1, pos, ref, alt);

        altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        impact = classifier.classifyVariant(var, negTrans);

        assertEquals(16, impact.codingContext().CodingBase);
        assertEquals(-3, impact.codingContext().NearestExonDistance);

        assertEquals("c.16-4_16-3dupAA", impact.codingContext().Hgvs);
    }
}
