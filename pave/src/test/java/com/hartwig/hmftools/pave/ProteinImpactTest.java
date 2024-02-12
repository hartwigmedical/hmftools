package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.codon.Codons.START_AMINO_ACID;
import static com.hartwig.hmftools.common.codon.Codons.START_CODON;
import static com.hartwig.hmftools.common.codon.Codons.STOP_CODON_1;
import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.START_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.STOP_GAINED;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.STOP_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SYNONYMOUS;
import static com.hartwig.hmftools.pave.ImpactTestUtils.generateTestBases;
import static com.hartwig.hmftools.pave.ImpactTestUtils.getAminoAcidCodon;
import static com.hartwig.hmftools.pave.ImpactTestUtils.getAminoAcidsCodons;
import static com.hartwig.hmftools.pave.impact.ProteinUtils.getExtraBases;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.fusion.FusionCommon;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.pave.impact.ImpactClassifier;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

import org.junit.Test;

public class ProteinImpactTest
{
    private final MockRefGenome mRefGenome;
    private final String mRefBases;
    private final ImpactClassifier mClassifier;
    private final TranscriptData mPosTrans;
    private final TranscriptData mNegTrans;

    public ProteinImpactTest()
    {
        mRefGenome = new MockRefGenome();

        // construct an exon with specific amino acids
        int preGene = 10;
        int prePostCoding = 10;
        String refBases = generateTestBases(preGene);

        // pos codons: M 20-22, A 23-25, C 26-28, D 29-31, L 32-34, L 35-37, G 38-40, H 41-43, E 44-46, stopX 47-49
        // M  A  C  D  L  L  G  H  E  X
        // ATGGCTTGTGACTTATTAGGACACGAGTAA

        refBases += generateTestBases(prePostCoding);
        String codingBases = getAminoAcidCodon(START_AMINO_ACID);
        codingBases += getAminoAcidsCodons("ACDLLGHE", false);
        codingBases += STOP_CODON_1;

        refBases += codingBases;
        refBases += generateTestBases(preGene);

        int transStart = preGene;
        int codingStart = transStart + prePostCoding;
        int codingEnd = codingStart + codingBases.length() - 1;
        int transEnd = codingEnd + prePostCoding;

        mPosTrans = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] {transStart}, transEnd - transStart,
                codingStart, codingEnd, false, "");

        int nextRand = 100 - transEnd;
        refBases += generateTestBases(nextRand - 1);

        refBases += generateTestBases(preGene);
        refBases += reverseComplementBases(codingBases);
        refBases += generateTestBases(preGene);

        transStart = transEnd + nextRand;
        codingStart = transStart + prePostCoding;
        codingEnd = codingStart + codingBases.length() - 1;
        transEnd = codingEnd + prePostCoding;

        // neg codons: M 139-37, A 136-134, C 133-131, D 130-128, L 127-125, L 124-122, G 121-119, H 118-116, E 115-113, stopX 112-110
        //   X  E  H  G  L  L  D  C  A  M
        // TTACTCGTGTCCTAATAAGTCACAAGCCAT GATCGATCGA
        // 110                            140
        mNegTrans =  createTransExons(
                GENE_ID_2, TRANS_ID_2, NEG_STRAND, new int[] {transStart}, transEnd - transStart,
                codingStart, codingEnd, false, "");

        mRefBases = refBases;
        mRefGenome.RefGenomeMap.put(CHR_1, refBases);

        mClassifier = new ImpactClassifier(mRefGenome);
    }

    @Test
    public void testSynonymousMissense()
    {
        // SNV - synonymous on first A, 3rd codon base
        int pos = 25;
        String ref = mRefBases.substring(pos, pos + 1);
        String alt = "C";
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = mClassifier.classifyVariant(var, mPosTrans);
        assertEquals(SYNONYMOUS, impact.topEffect());

        assertEquals("c.6T>C", impact.codingContext().Hgvs);
        assertEquals("p.Ala2=", impact.proteinContext().Hgvs);

        checkHgvsStrings(pos, 1, alt, SYNONYMOUS, "c.6T>C", "p.Ala2=");

        // same codon but missense
        pos = 24;
        alt = "G";
        checkHgvsStrings(pos, 1, alt, MISSENSE, "c.5C>G", "p.Ala2Gly");

        // MNV same codon synonymous: L=TTA -> CTG
        pos = 32;
        alt = "CTG";

        checkHgvsStrings(pos, 3, alt, SYNONYMOUS, "c.13_15delTTAinsCTG", "p.Leu5=");

        // now missense
        alt = "GGC";
        checkHgvsStrings(pos, 3, alt, MISSENSE, "c.13_15delTTAinsGGC", "p.Leu5Gly");

        // MNV synonymous spanning 2 codons: L L, TTA TTA -> TTG -> CTG
        pos = 34;
        alt = "GC";
        checkHgvsStrings(pos, 2, alt, SYNONYMOUS, "c.15_16delATinsGC", "p.Leu5_Leu6=");

        // now missense: L L, TTA TTA -> TCG -> GTA
        pos = 33;
        alt = "CGG";
        checkHgvsStrings(pos, 3, alt, MISSENSE, "c.14_16delTATinsCGG", "p.Leu5_Leu6delinsSerVal");

        // causing stop gained, across 2 codons, L(TTA) L(TTA) -> (C)TGT (X)TAA
        // eg p.Cys495_Val496delinsArg*
        pos = 33;
        alt = "GTTA";
        checkHgvsStrings(pos, 4, alt, STOP_GAINED, "c.14_17delTATTinsGTTA", "p.Leu5_Leu6delinsCys*");

        // start lost
        pos = 21;
        alt = "GT";
        checkHgvsStrings(pos, 2, alt, START_LOST, "c.2_3delTGinsGT", "p.Met1?");

        // stop lost
        pos = 48;
        alt = "C";
        checkHgvsStrings(pos, 1, alt, STOP_LOST, "c.29A>C", "p.Ter10Serext*?");

        // stop gained but in the first codon where an MNV is affecting 2
        pos = 27;
        alt = "AAC";
        checkHgvsStrings(pos, 3, alt, STOP_GAINED, "c.8_10delGTGinsAAC", "p.Cys3_Asp4delins*");
    }

    @Test
    public void testInframeDeletion()
    {
        // p.Lys2del
        // conservative DEL: codon 2 = A
        int pos = 22;
        String ref = mRefBases.substring(pos, pos + 4);
        String alt = ref.substring(0, 1);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = mClassifier.classifyVariant(var, mPosTrans);
        assertEquals(INFRAME_DELETION, impact.topEffect());

        assertEquals("c.4_6delGCT", impact.codingContext().Hgvs);
        assertEquals("p.Ala2del", impact.proteinContext().Hgvs);

        checkHgvsStrings(pos, 4, alt, INFRAME_DELETION, "c.4_6delGCT", "p.Ala2del");

        // 3 codons conservative
        checkHgvsStrings(pos, 10, alt, INFRAME_DELETION, "c.4_12delGCTTGTGAC", "p.Ala2_Asp4del");

        // 2 codons non-conservative: deletes part of A and C and makes a G
        pos = 23;
        alt = mRefBases.substring(pos, pos + 1);

        checkHgvsStrings(pos, 4, alt, INFRAME_DELETION, "c.5_7delCTT", "p.Ala2_Cys3delinsGly");

        // causing a stop gained: L(TTA) G(GGA) -> TGA
        pos = 35;
        alt = mRefBases.substring(pos, pos + 1);
        checkHgvsStrings(pos, 4, alt, STOP_GAINED, "c.17_19delTAG", "p.Leu6_Gly7delins*");

        // causing a stop lost H(CAC) E(GAG) X(TAA) but report as frameshift since the first AA to change is not a stop
        pos = 41;
        alt = mRefBases.substring(pos, pos + 1);
        checkHgvsStrings(pos, 7, alt, STOP_LOST, "c.23_28delACGAGT", "p.His8fs");

        // causing a stop lost E(GAG) X(TAA) but the inframe DEL makes a synonymous E
        pos = 44;
        alt = mRefBases.substring(pos, pos + 1);
        checkHgvsStrings(pos, 4, alt, STOP_LOST, "c.26_28delAGT", "p.Glu9ext*?");
    }

    @Test
    public void testInframeInsertion()
    {
        // duplication (single AA) p.Gln8dup
        // duplication (range) p.Gly4_Gln6dup
        // insertion p.Lys2_Leu3insGlnSer
        // insertion (conservative stop) p.Ser81_Val82ins*
        // insertion (non conservative) p.Cys28delinsTrpVal

        // insert a codon in between A1 and C2
        int pos = 25;
        String ref = mRefBases.substring(pos, pos + 1);
        String alt = ref + getAminoAcidCodon('R');
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = mClassifier.classifyVariant(var, mPosTrans);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        assertEquals("c.6_7insAGG", impact.codingContext().Hgvs);
        assertEquals("p.Ala2_Cys3insArg", impact.proteinContext().Hgvs);

        checkHgvsStrings(pos, 1, alt, INFRAME_INSERTION, "c.6_7insAGG", "p.Ala2_Cys3insArg");

        // insert multiple codons
        alt = ref + getAminoAcidsCodons("RVS", false);
        checkHgvsStrings(pos, 1, alt, INFRAME_INSERTION, "c.6_7insAGGGTGTCA", "p.Ala2_Cys3insArgValSer");

        // non-conservative insert: CDL -> CEGCL   C + D(GAC) + L(TTA) + L -> C + E(GAA) (G)GGG (C)TGT + L
        pos = 30;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + getAminoAcidsCodons("RV", false);
        checkHgvsStrings(pos, 1, alt, INFRAME_INSERTION, "c.11_12insAGGGTG", "p.Asp4delinsGluGlyCys");

        // pos codons: M 20-22, A 23-25, C 26-28, D 29-31, L 32-34, L 35-37, G 38-40, H 41-43, E 44-46, stopX 47-49
        // M  A  C  D  L  L  G  H  E  X
        // ATGGCTTGTGACTTATTAGGACACGAGTAA

        // DLLG -> DCCLLG, or shortened DL -> DCCL, so should be p.D4_L5insCC
        // ref: L
        //      T       TA
        // alt: C     C   L
        //      T  GT TGT T  TA
        pos = 32;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + "GTTGTT"; // getAminoAcidsCodons("CC", false);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = mClassifier.classifyVariant(var, mPosTrans);
        assertEquals(INFRAME_INSERTION, impact.topEffect());
        assertEquals("p.Asp4_Leu5insCysCys", impact.proteinContext().Hgvs);

        // neg codons: M 139-37, A 136-134, C 133-131, D 130-128, L 127-125, L 124-122, G 121-119, H 118-116, E 115-113, stopX 112-110
        //   X  E  H  G  L  L  D  C  A  M
        // TTACTCGTGTCCTAATAAGTCACAAGCCAT GATCGATCGA
        // 110                            140

        pos = 123;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + getAminoAcidsCodons("GC", true);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = mClassifier.classifyVariant(var, mNegTrans);
        assertEquals(INFRAME_INSERTION, impact.topEffect());
        // assertEquals("p.Leu5_Leu6insCysGly", impact.proteinContext().Hgvs);

        // duplications

        // single G repeated
        pos = 37;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + "GGA";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = mClassifier.classifyVariant(var, mPosTrans);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        assertEquals("c.18_19insGGA", impact.codingContext().Hgvs);
        assertEquals("p.Gly7dup", impact.proteinContext().Hgvs);

        // repeated more than 1 is no longer a dup
        pos = 37;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + "GGAGGAGGA";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = mClassifier.classifyVariant(var, mPosTrans);
        assertEquals(INFRAME_INSERTION, impact.topEffect());
        assertEquals("p.Leu6_Gly7insGlyGlyGly", impact.proteinContext().Hgvs);

        // range of AAs: L(TTA) G(GGA) H(CAC) E(GAG), T>TAGGACACGA - duplicating GHE
        pos = 36;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + "AGGACACGA";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = mClassifier.classifyVariant(var, mPosTrans);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        assertEquals("c.17_18insAGGACACGA", impact.codingContext().Hgvs);
        assertEquals("p.Gly7_Glu9dup", impact.proteinContext().Hgvs);

        // both Ls are duplicated but treat this is a range rather than repeat of a single
        pos = 31;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + "TTATTA";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = mClassifier.classifyVariant(var, mPosTrans);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        assertEquals("p.Leu5_Leu6dup", impact.proteinContext().Hgvs);
    }

    @Test
    public void testFrameshifts()
    {
        int pos = 26;
        String ref = mRefBases.substring(pos, pos + 2);
        String alt = ref.substring(0, 1);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = mClassifier.classifyVariant(var, mPosTrans);
        assertEquals(FRAMESHIFT, impact.topEffect());

        assertEquals("c.8delG", impact.codingContext().Hgvs);
        assertEquals("p.Cys3fs", impact.proteinContext().Hgvs);

        checkHgvsStrings(pos, 2, alt, FRAMESHIFT, "c.8delG", "p.Cys3fs");

        // longer DEL still uses first ref codon
        checkHgvsStrings(pos, 12, alt, FRAMESHIFT, "c.8_18delGTGACTTATTA", "p.Cys3fs");

        // from an insert
        alt = mRefBases.substring(pos, pos + 1) + "CC";
        checkHgvsStrings(pos, 1, alt, FRAMESHIFT, "c.7_8insCC", "p.Cys3fs");

        // causing a stop gained: L(TTA) -> TGA, and no stop-gained identifier even if one is added
        pos = 35;
        ref = mRefBases.substring(pos, pos + 1);
        alt = ref + "GA";
        checkHgvsStrings(pos, 1, alt, FRAMESHIFT, "c.16_17insGA", "p.Leu6*");

        // SNV and deletion combined: ie CT -> A, should be c.5_6delCTinsA instead of c.6delT
        pos = 24;
        ref = mRefBases.substring(pos, pos + 2);
        alt = "A";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = mClassifier.classifyVariant(var, mPosTrans);
        assertEquals(FRAMESHIFT, impact.topEffect());

        assertEquals(5, impact.codingContext().CodingBase);
        assertEquals(1, impact.codingContext().DeletedCodingBases);
        assertEquals("c.5_6delCTinsA", impact.codingContext().Hgvs);
        assertEquals("p.Ala2fs", impact.proteinContext().Hgvs);
    }

    @Test
    public void testSpecificScenarios()
    {
        int preGene = 10;
        int prePostCoding = 10;
        String refBases = generateTestBases(preGene);

        // codons: M 20-22, G 23-25, P 26-28, P 29-31, H 32-34, L 35-37, stopX
        // M  G  P  H  L  X
        // ATGGGACCCCACTTATAA
        // 20

        refBases += generateTestBases(prePostCoding);
        String codingBases = getAminoAcidCodon(START_AMINO_ACID);
        codingBases += "GGACCCCCCCACTTA"; // G, P, P, H, L
        codingBases += STOP_CODON_1;

        refBases += codingBases;
        refBases += generateTestBases(preGene);

        mRefGenome.RefGenomeMap.put(CHR_2, refBases);

        int transStart = preGene;
        int codingStart = transStart + prePostCoding;
        int codingEnd = codingStart + codingBases.length() - 1;
        int transEnd = codingEnd + prePostCoding;

        TranscriptData transData = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] {transStart}, transEnd - transStart,
                codingStart, codingEnd, false, "");

        // frameshift DEL
        int pos = 25;
        String ref = refBases.substring(pos, pos + 2);
        String alt = ref.substring(0, 1);
        VariantData var = new VariantData(CHR_2, pos, ref, alt);

        VariantTransImpact impact = mClassifier.classifyVariant(var, transData);
        assertEquals(FRAMESHIFT, impact.topEffect());

        assertEquals("c.7delC", impact.codingContext().Hgvs);
        assertEquals("p.His5fs", impact.proteinContext().Hgvs);

        // conservative inframe deletion with homology - deletes a P, but can be shifted further
        pos = 25;
        ref = refBases.substring(pos, pos + 4);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_2, pos, ref, alt);

        impact = mClassifier.classifyVariant(var, transData);
        assertEquals(INFRAME_DELETION, impact.topEffect());

        assertEquals("c.7_9delCCC", impact.codingContext().Hgvs);
        assertEquals("p.Pro4del", impact.proteinContext().Hgvs);
    }

    @Test
    public void testSpecificScenarios2()
    {
        int preGene = 10;
        int prePostCoding = 10;
        String refBases = generateTestBases(preGene);

        // codons: X 20-22, K 23-25, K 26-28, K 29-31, M 32-34
        // ATT TTTTTTTGC CAT
        // 20

        refBases += generateTestBases(prePostCoding);
        String codingBases = "ATTTTTTTTTGCCAT"; // M, A, K, K, X

        refBases += codingBases;
        refBases += generateTestBases(preGene);

        mRefGenome.RefGenomeMap.put(CHR_2, refBases);

        int transStart = preGene;
        int codingStart = transStart + prePostCoding;
        int codingEnd = codingStart + codingBases.length() - 1;
        int transEnd = codingEnd + prePostCoding;

        TranscriptData transData = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, new int[] {transStart}, transEnd - transStart,
                codingStart, codingEnd, false, "");

        int pos = 28;
        String ref = refBases.substring(pos, pos + 1);
        String alt = ref.substring(0, 1) + "T";
        VariantData var = new VariantData(CHR_2, pos, ref, alt);

        VariantTransImpact impact = mClassifier.classifyVariant(var, transData);
        assertEquals(FRAMESHIFT, impact.topEffect());

        assertEquals("c.6_7insA", impact.codingContext().Hgvs);
        // assertEquals("p.Ala5fs", impact.proteinContext().Hgvs);

        // test 2: inframe deletion of 2 codons leading up to exon boundary
        refBases = generateTestBases(preGene);

        refBases += generateTestBases(preGene);
        refBases += reverseComplementBases(STOP_CODON_1);
        refBases += getAminoAcidsCodons("GHC", true);
        refBases += generateTestBases(preGene); // intron
        refBases += getAminoAcidsCodons("LRD", true);
        refBases += reverseComplementBases(START_CODON);
        refBases += generateTestBases(preGene);

        int codingLength = 8 * 3;

        transStart = preGene;
        codingStart = transStart + prePostCoding;
        codingEnd = 53;

        // 10-19 3'UTR, 20-22 X, 23-25 C, 26-28 H 29-31 G, 32-41 intron 42-44 D, 45-47 R, 48-50 L, 51-53 M, 54-63 5'UTR
        // TTACTCGTGTCCTAATAAGTCACAAGCCAT GATCGATCGA
        // 110                            140
        int[] exonStarts = new int[] {10, 42};

        TranscriptData negTrans =  createTransExons(
                GENE_ID_2, TRANS_ID_2, NEG_STRAND, exonStarts, 21,
                codingStart, codingEnd, false, "");

        mRefGenome.RefGenomeMap.put(CHR_3, refBases);

        // deleting last 2 amino acids of second coding exon
        pos = 25;
        ref = refBases.substring(pos, pos + 7);
        alt = ref.substring(0, 1);
        var = new VariantData(CHR_3, pos, ref, alt);

        impact = mClassifier.classifyVariant(var, negTrans);
        assertEquals(INFRAME_DELETION, impact.topEffect());
        assertEquals("c.13_18delGGACAC", impact.codingContext().Hgvs);
        assertEquals("p.Gly5_His6del", impact.proteinContext().Hgvs);
    }

    private void checkHgvsStrings(
        int pos, int refLen, final String alt, final VariantEffect effect, final String codingStr, final String proteinStr)
    {
        // check positive then negative strand
        String ref = mRefBases.substring(pos, pos + refLen);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = mClassifier.classifyVariant(var, mPosTrans);
        assertEquals(effect, impact.topEffect());

        assertEquals(codingStr, impact.codingContext().Hgvs);
        assertEquals(proteinStr, impact.proteinContext().Hgvs);

        int negPosEquiv = mNegTrans.CodingEnd - (pos - mPosTrans.CodingStart);
        int negPos = refLen == alt.length() ? negPosEquiv - (refLen - 1) : negPosEquiv - refLen;
        String negRef = mRefBases.substring(negPos, negPos + refLen);

        String negAlt;

        if(var.isInsert())
            negAlt = negRef + reverseComplementBases(alt.substring(1));
        else if(var.isDeletion())
            negAlt = negRef.substring(0, 1);
        else
            negAlt = reverseComplementBases(alt);

        VariantData varNeg = new VariantData(CHR_1, negPos, negRef, negAlt);

        VariantTransImpact impactNeg = mClassifier.classifyVariant(varNeg, mNegTrans);
        assertEquals(effect, impactNeg.topEffect());

        assertEquals(codingStr, impactNeg.codingContext().Hgvs);
        assertEquals(proteinStr, impactNeg.proteinContext().Hgvs);
    }

    @Test
    public void testSynonymousMissenseImpactsOld()
    {
        // inferior version of the above, remove once new tests are complete
        final MockRefGenome refGenome = new MockRefGenome();

        int[] exonStarts = { 0, 100, 200 };

        // codons start on at 10, 13, 16 etc
        Integer codingStart = Integer.valueOf(10);
        Integer codingEnd = Integer.valueOf(250);

        TranscriptData transDataPos = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 80, codingStart, codingEnd, false, "");

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        String chr1Bases = generateTestBases(300);

        // set the specific AAs for a region of this mock ref genome
        // S: TCA -> TCG - so last base of codon can change
        // I: ATC -> ATT
        // L: TTA -> CTA - first base changes

        //                     40      43      46
        String aminoAcidSeq = "TCA" + "ATC" + "TTA";
        chr1Bases = chr1Bases.substring(0, 40) + aminoAcidSeq + chr1Bases.substring(49);
        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        // test SNVs at the 1st and 3rd codon bases

        // last base of a codon changes
        int codonPos = 40;
        String codon = chr1Bases.substring(codonPos, codonPos + 3);
        String aminoAcid = AminoAcids.findAminoAcidForCodon(codon);

        String alt = "G";
        String synCodon = chr1Bases.substring(codonPos, codonPos + 2) + alt;
        assertTrue(aminoAcid.equals(AminoAcids.findAminoAcidForCodon(synCodon)));

        int pos = 42;
        String ref = chr1Bases.substring(pos, pos + 1);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        VariantTransImpact impact = classifier.classifyVariant(var, transDataPos);
        assertEquals(SYNONYMOUS, impact.topEffect());

        // now missense
        pos = 41;
        ref = chr1Bases.substring(pos, pos + 1);
        alt = "T";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);
        assertEquals(MISSENSE, impact.topEffect());

        // first base of a codon changes
        codonPos = 46;
        codon = chr1Bases.substring(codonPos, codonPos + 3);
        aminoAcid = AminoAcids.findAminoAcidForCodon(codon);

        alt = "C";
        synCodon = alt + chr1Bases.substring(codonPos + 1, codonPos + 3);
        assertTrue(aminoAcid.equals(AminoAcids.findAminoAcidForCodon(synCodon)));

        pos = 46;
        ref = chr1Bases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);
        assertEquals(SYNONYMOUS, impact.topEffect());

        // test for an MNV spanning 3 codons - first in the last codon pos at 42 then all the next and 2 into the final one at 46
        pos = 42;
        ref = chr1Bases.substring(pos, pos + 5);
        alt = "G" + "ATT" + "C";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);
        assertEquals(SYNONYMOUS, impact.topEffect());

        // now missense by change middle codon
        alt = "G" + "AAT" + "C";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataPos);
        assertEquals(MISSENSE, impact.topEffect());

        // test reverse strand
        TranscriptData transDataNeg = createTransExons(
                GENE_ID_2, TRANS_ID_2, FusionCommon.NEG_STRAND, exonStarts, 80, codingStart, codingEnd, false, "");

        // coding starts at 250 so codons start at 250, 247, 244 etc
        // still change the AA sequence S, I then L - for the range 239-241, 242-244 and 245-247
        String aminoAcidSeqRev = reverseComplementBases(aminoAcidSeq);
        chr1Bases = chr1Bases.substring(0, 239) + aminoAcidSeqRev + chr1Bases.substring(248);
        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        codonPos = 242;
        codon = chr1Bases.substring(codonPos, codonPos + 3);
        aminoAcid = AminoAcids.findAminoAcidForCodon(reverseComplementBases(codon));

        // change last base of a codon, which is the first base
        alt = "A";
        synCodon = alt + chr1Bases.substring(codonPos + 1, codonPos + 3);
        assertTrue(aminoAcid.equals(AminoAcids.findAminoAcidForCodon(reverseComplementBases(synCodon))));

        pos = 242;
        ref = chr1Bases.substring(pos, pos + 1);
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);
        assertEquals(SYNONYMOUS, impact.topEffect());

        // now missense
        pos = 243;
        ref = chr1Bases.substring(pos, pos + 1);
        alt = "T";
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);
        assertEquals(MISSENSE, impact.topEffect());

        // test a MNV spanning 3 codons as before
        pos = 241;
        ref = chr1Bases.substring(pos, pos + 5);
        alt = reverseComplementBases("G" + "ATT" + "C");
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);
        assertEquals(SYNONYMOUS, impact.topEffect());

        // and missense
        alt = reverseComplementBases("G" + "AAT" + "C");
        var = new VariantData(CHR_1, pos, ref, alt);

        impact = classifier.classifyVariant(var, transDataNeg);
        assertEquals(MISSENSE, impact.topEffect());
    }

    @Test
    public void testExtraCodingBases()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        String chr1Bases = generateTestBases(100);
        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        int[] exonStarts = { 10, 30, 50, 70};

        // codons start on at 10, 13, 16 etc
        Integer codingStart = Integer.valueOf(15);
        Integer codingEnd = Integer.valueOf(75);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        ExonData exon = transDataPosStrand.exons().get(1);

        String bases = getExtraBases(transDataPosStrand, refGenome, CHR_1, exon, 32, 1, false);
        String refBases = chr1Bases.substring(31, 32);
        assertTrue(refBases.equals(bases));

        bases = getExtraBases(transDataPosStrand, refGenome, CHR_1, exon, 32, 2, false);
        refBases = chr1Bases.substring(30, 32);
        assertTrue(refBases.equals(bases));

        bases = getExtraBases(transDataPosStrand, refGenome, CHR_1, exon, 32, 3, false);
        refBases = chr1Bases.substring(20, 21) + chr1Bases.substring(30, 32);
        assertTrue(refBases.equals(bases));

        bases = getExtraBases(transDataPosStrand, refGenome, CHR_1, exon, 38, 4, true);
        refBases = chr1Bases.substring(39, 41) + chr1Bases.substring(50, 52);
        assertTrue(refBases.equals(bases));

        TranscriptData transDataNegStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        exon = transDataPosStrand.exons().get(2);

        bases = getExtraBases(transDataPosStrand, refGenome, CHR_1, exon, 52, 4, false);
        refBases = chr1Bases.substring(39, 41) + chr1Bases.substring(50, 52);
        assertTrue(refBases.equals(bases));

        bases = getExtraBases(transDataPosStrand, refGenome, CHR_1, exon, 58, 6, true);
        refBases = chr1Bases.substring(59, 61) + chr1Bases.substring(70, 74);
        assertTrue(refBases.equals(bases));

    }

}
