package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_INSERTION;
import static com.hartwig.hmftools.pave.ImpactTestUtils.generateTestBases;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;
import static com.hartwig.hmftools.pave.impact.PaveUtils.createRightAlignedVariant;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.FusionCommon;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.pave.impact.CodingContext;
import com.hartwig.hmftools.pave.impact.HgvsCoding;
import com.hartwig.hmftools.pave.impact.ImpactClassifier;
import com.hartwig.hmftools.pave.impact.PaveUtils;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

import org.junit.Ignore;
import org.junit.Test;

public class HgvsCodingTest
{
    private final GeneDataCache GENE_DATA_CACHE = new GeneDataCache("", V37, null);
    private final MockRefGenome REF_GENOME = new MockRefGenome();
    private final ImpactClassifier CLASSIFIER = new ImpactClassifier(REF_GENOME);

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
                GENE_ID_1, TRANS_ID_1, FusionCommon.POS_STRAND, new int[] { 10, 50 }, 24,
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
                GENE_ID_1, TRANS_ID_1, FusionCommon.NEG_STRAND, new int[] { 10, 50 }, 24,
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
        assertEquals(-2, impact.codingContext().NearestExonDistance); // nearest exon ends at 34, using upstream pos of 37

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

    private void addGeneTranscript(final TranscriptData transcriptData)
    {
        GeneData geneData = createEnsemblGeneData(
                GENE_ID_1, GENE_NAME_1, CHR_1, transcriptData.Strand, transcriptData.TransStart, transcriptData.TransEnd);

        // clear for each test
        GENE_DATA_CACHE.getEnsemblCache().getChrGeneDataMap().clear();
        GENE_DATA_CACHE.getEnsemblCache().getTranscriptDataMap().clear();

        GENE_DATA_CACHE.getEnsemblCache().getChrGeneDataMap().put(geneData.Chromosome, Lists.newArrayList(geneData));
        GENE_DATA_CACHE.getEnsemblCache().getTranscriptDataMap().put(geneData.GeneId, Lists.newArrayList(transcriptData));
    }

    @Test
    public void testAminoAcidRightAlignment()
    {
        // while there is no homology, the amino acid impact is right aligned - only relevant for positive-strand genes

        // amino acids:           M   A   P   P   P   P   Q
        // AA index:              1   2   3   4   5   6   7
        // exon 1:                20           30
        // position:              012 345 678 901 234 567 890
        // coding index:          123 456 789 012 345 678 901
        String refCodingBases1 = "ATG GCT CCA CCC CCG CCT CAG".replace(" ", "");
        String intronicBases = "TTTTGGGGCCCCAAA"; // as in previous test

        String refBases = generateTestBases(20) + refCodingBases1 + intronicBases;
        REF_GENOME.RefGenomeMap.put(CHR_1, refBases);

        TranscriptData posTrans = createTransExons(
                GENE_ID_1, TRANS_ID_1, FusionCommon.POS_STRAND, new int[] { 10, 56 }, 30,
                20, 70, false, "");

        addGeneTranscript(posTrans);

        // duplication of first P should be converted to a right-aligned AA change
        int pos = 25;
        String ref = refBases.substring(pos, pos + 1);
        String alt = refBases.substring(pos, pos + 4);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        String altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        PaveUtils.findVariantImpacts(var, CLASSIFIER, GENE_DATA_CACHE, null);

        assertTrue(var.getImpacts().containsKey(GENE_NAME_1));
        VariantTransImpact impact = var.getImpacts().get(GENE_NAME_1).get(0);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        assertEquals("c.7_9dupCCA", impact.codingContext().Hgvs);
        assertEquals("p.Pro6dup", impact.proteinContext().Hgvs);
    }

    @Test
    public void testAminoAcidRightAlignment2()
    {
        // amino acids:         M   Q   L   L   R   G
        // AA index:            1   2   3   4   5   6
        // exon 1:              20           30
        // position:            012 345 678 901
        // coding index:        123 456 789 012 345 678
        String refCodingBases ="ATG CAA TTA TTA AGG GGA CAC".replace(" ", "");
        String intronicBases = "TTTTGGGGCCCCAAA";

        // variant is inserting RV = AGA (R) and GTG (V)
        // ref: A         AGG
        // alt: A AGA GTG AGG

        MockRefGenome refGenome = new MockRefGenome();

        String refBases = generateTestBases(20) + refCodingBases + intronicBases;
        refGenome.RefGenomeMap.put(CHR_1, refBases);

        GeneDataCache geneDataCache = new GeneDataCache("", V37, null);

        GeneData geneData = createEnsemblGeneData(GENE_ID_1, GENE_NAME_1, CHR_1, FusionCommon.POS_STRAND, 20, 90);

        geneDataCache.getEnsemblCache().getChrGeneDataMap().put(CHR_1, Lists.newArrayList(geneData));

        TranscriptData posTrans = createTransExons(
                GENE_ID_1, TRANS_ID_1, FusionCommon.POS_STRAND, new int[] { 10, 56 }, 30,
                20, 70, false, "");

        geneDataCache.getEnsemblCache().getTranscriptDataMap().put(GENE_ID_1, Lists.newArrayList(posTrans));

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        // duplication of first P should be converted to a right-aligned AA change
        int pos = 31;
        String ref = refBases.substring(pos, pos + 1);
        String alt = ref + "AGAGTG";
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        String altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        PaveUtils.findVariantImpacts(var, classifier, geneDataCache, null);

        assertTrue(var.getImpacts().containsKey(GENE_NAME_1));
        VariantTransImpact impact = var.getImpacts().get(GENE_NAME_1).get(0);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        // assertEquals("c.12insAGAGTG", impact.codingContext().Hgvs);
        assertEquals("p.Arg5_Gly6insValArg", impact.proteinContext().Hgvs);
    }

    @Test
    public void testAminoAcidRightAlignment3()
    {
        // amino acids:           M   G   L   *
        // AA index:              1   2   3   4
        // exon 1:                20           30
        // position:              012 345 678 901
        // coding index:          123 456 789 012
        String refCodingBases =  "ATG GGA CTA TGA".replace(" ", "");
        String intronicBases = "TTTTGGGGCCCCAAA";

        MockRefGenome refGenome = new MockRefGenome();

        String refBases = generateTestBases(20) + refCodingBases + intronicBases;
        refGenome.RefGenomeMap.put(CHR_1, refBases);

        GeneDataCache geneDataCache = new GeneDataCache("", V37, null);

        GeneData geneData = createEnsemblGeneData(GENE_ID_1, GENE_NAME_1, CHR_1, FusionCommon.POS_STRAND, 20, 90);

        geneDataCache.getEnsemblCache().getChrGeneDataMap().put(CHR_1, Lists.newArrayList(geneData));

        TranscriptData posTrans = createTransExons(
                GENE_ID_1, TRANS_ID_1, FusionCommon.POS_STRAND, new int[] { 10, 56 }, 30,
                20, 70, false, "");

        geneDataCache.getEnsemblCache().getTranscriptDataMap().put(GENE_ID_1, Lists.newArrayList(posTrans));

        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        // duplication of first P should be converted to a right-aligned AA change
        int pos = 25;
        String ref = refBases.substring(pos, pos + 2);
        String alt = refBases.substring(pos, pos + 1);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        String altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        PaveUtils.findVariantImpacts(var, classifier, geneDataCache, null);

        assertTrue(var.getImpacts().containsKey(GENE_NAME_1));
        VariantTransImpact impact = var.getImpacts().get(GENE_NAME_1).get(0);
        assertEquals(FRAMESHIFT, impact.topEffect());

        assertEquals("c.7delC", impact.codingContext().Hgvs);
        assertEquals("p.Leu3fs", impact.proteinContext().Hgvs);
    }

    @Ignore
    @Test
    public void testRightAlignmentOfCodingOfInFrameDuplication()
    {
        // pos codons: M 20-22, A 23-25, D 26-28, A 29-31, D 32-34, S 35-37, Q 38-40, L 41-43, G 44-46, H 47-49, E 50-52, stopX 47-49
        // amino acids:           M   A   D   A   D   S   Q
        // exon 1:                20           30
        // position:              012 345 678 901 234 567 890
        String refCodingBases1 = "ATG GCT GAT GCT GAT TCG CAG".replace(" ", "");

        // intron:              41     50
        // position             123456789012345
        // bases into intron:   12345   7654321
        String intronicBases = "TTTTGGGGCCCCAAA";

        // amino acids:           L   G   H   E   X
        // exon 2:                56   60           70
        // position:              678 901 234 567 890
        String refCodingBases2 = "TTA GGA CAC GAG TAA".replace(" ", "");

        MockRefGenome refGenome = new MockRefGenome();

        String refBases = generateTestBases(20) + refCodingBases1 + intronicBases + refCodingBases2 + generateTestBases(20);
        refGenome.RefGenomeMap.put(CHR_1, refBases);
        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        TranscriptData posTrans = createTransExons(
                GENE_ID_1, TRANS_ID_1, FusionCommon.POS_STRAND, new int[] {10, 56}, 30,
                20, 70, false, "");

        addGeneTranscript(posTrans);

        // duplication of first AD in ADAD should be reported right-maximal
        // ATG [GCT GAT] GCT GAT TCG CAG -> ATG [GCT GAT GCT GAT] GCT GAT TCG CAG = ATG GCT GAT [GCT GAT GCT GAT] TCG CAG
        int pos = 23;
        String ref = refBases.substring(pos, pos + 1);
        String alt = refBases.substring(pos, pos + 7);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        String altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        var.setRealignedVariant(createRightAlignedVariant(var, REF_GENOME));

        // TODO: check if use of RA var passes
        PaveUtils.findVariantImpacts(var, CLASSIFIER, GENE_DATA_CACHE, null);

        VariantTransImpact impact = var.getImpacts().get(GENE_NAME_1).get(0);
        // VariantTransImpact impact = classifier.classifyVariant(var, posTrans);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        assertEquals("c.9_14dupGCTGAT", impact.codingContext().Hgvs);
        assertEquals("p.Ala4_Asp5dup", impact.proteinContext().Hgvs);

        // duplication of DA in ADAD should give the same result
        // ATG GCT [GAT GCT] GAT TCG CAG -> ATG GCT [GAT GCT GAT GCT] GAT TCG CAG = ATG GCT GAT [GCT GAT GCT GAT] TCG CAG
        pos = 23 + 3;
        ref = refBases.substring(pos, pos + 1);
        alt = refBases.substring(pos, pos + 7);
        var = new VariantData(CHR_1, pos, ref, alt);

        altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        PaveUtils.findVariantImpacts(var, CLASSIFIER, GENE_DATA_CACHE, null);
        impact = var.getImpacts().get(GENE_NAME_1).get(0);

        // impact = classifier.classifyVariant(var, posTrans);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        assertEquals("c.9_14dupGCTGAT", impact.codingContext().Hgvs);
        assertEquals("p.Ala4_Asp5dup", impact.proteinContext().Hgvs);

        // duplication of second AD in ADAD should give the same result
        // ATG GCT GAT [GCT GAT] TCG CAG -> ATG GCT GAT [GCT GAT GCT GAT] TCG CAG
        pos = 23 + 6;
        ref = refBases.substring(pos, pos + 1);
        alt = refBases.substring(pos, pos + 7);
        var = new VariantData(CHR_1, pos, ref, alt);

        altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        impact = classifier.classifyVariant(var, posTrans);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        assertEquals("c.9_14dupGCTGAT", impact.codingContext().Hgvs);
        assertEquals("p.Ala4_Asp5dup", impact.proteinContext().Hgvs);
    }

    @Ignore
    @Test
    public void testRightAlignmentOfCodingOfFrameshiftDuplication()
    {
        // pos codons: M 20-22, A 23-25, D 26-28, M 29-31, A 32-34, S 35-37, Q 38-40, L 41-43, G 44-46, H 47-49, E 50-52, stopX 47-49
        // amino acids:           M   A   D   M   A   S   Q
        // exon 1:                20           30
        // position:              012 345 678 901 234 567 890
        String refCodingBases1 = "ATG GCT GAT ATG GCC TCG CAG".replace(" ", "");
        String intronicBases = "TTTTGGGGCCCCAAA"; // same as previous test
        String refCodingBases2 = "TTAGGACACGAGTAA"; // ditto

        MockRefGenome refGenome = new MockRefGenome();

        String refBases = generateTestBases(20) + refCodingBases1 + intronicBases + refCodingBases2 + generateTestBases(20);
        refGenome.RefGenomeMap.put(CHR_1, refBases);
        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        TranscriptData posTrans = createTransExons(
                GENE_ID_1, TRANS_ID_1, FusionCommon.POS_STRAND, new int[] { 10, 56 }, 30,
                20, 70, false, "");

        int pos = 29;
        String ref = refBases.substring(pos, pos + 1);
        String alt = refBases.substring(pos, pos + 6);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        String altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        VariantTransImpact impact = classifier.classifyVariant(var, posTrans);
        assertEquals(INFRAME_INSERTION, impact.topEffect());

        // ATG GCT GAT A[TG GCC] TCG CAG -> ATG GCT GAT A[TG GCC TG GCC] TCG CAG = ATG GCT GAT AT[G GCC TG GCC T]CG CAG
        assertEquals("c.9_14dupGCTGAT", impact.codingContext().Hgvs);
        // new AA seq is M A D M A W P
        assertEquals("p.Trp6fs", impact.proteinContext().Hgvs);
    }

    @Ignore
    @Test
    public void testRightAlignmentOfCodingOfFrameshiftDeletion()
    {
        // amino acids:           M   K   E   K   K   K   P
        // exon 1:                20           30
        // position:              012 345 678 901 234 567 890
        String refCodingBases1 = "ATG AAG GAA AAA AAA AAG CCT".replace(" ", "");
        String intronicBases = "TTTTGGGGCCCCAAA"; // same as previous test
        String refCodingBases2 = "TTAGGACACGAGTAA"; // ditto

        MockRefGenome refGenome = new MockRefGenome();

        String refBases = generateTestBases(20) + refCodingBases1 + intronicBases + refCodingBases2 + generateTestBases(20);
        refGenome.RefGenomeMap.put(CHR_1, refBases);
        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        TranscriptData posTrans = createTransExons(
                GENE_ID_1, TRANS_ID_1, FusionCommon.POS_STRAND, new int[] { 10, 56 }, 30,
                20, 70, false, "");

        int pos = 34;
        String ref = refBases.substring(pos, pos + 3);
        String alt = refBases.substring(pos, pos + 1);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        String altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);

        VariantTransImpact impact = classifier.classifyVariant(var, posTrans);
        assertEquals(FRAMESHIFT, impact.topEffect());

        assertEquals("c.16_17delAA", impact.codingContext().Hgvs);
        assertEquals("p.Lys6fs", impact.proteinContext().Hgvs);
    }

    @Test
    public void testDupIn5PUtr()
    {
        // random(20) exon1 intron1 exon2 intron2 exon3(coding) intron3 exon4(coding) random(20)
        // all exons and intronns of length 12
        String exon1 = "CTAGGACACGAG";
        String intron1 = "CGAGGGCCCAAA";
        String exon2 = "GGACACGAGTAA";
        String intron2 = "GGGCCCAAATTT";
        String exon3 = "ATGAAGGAACCT";
        String intron3 = "AAAGGGCCCTTT";
        String exon4 = "AAGATGGAACCT";
        String refBases = generateTestBases(20) + exon1 + intron1 + exon2 + intron2 + exon3 + intron3 + exon4 + generateTestBases(20);

        MockRefGenome refGenome = new MockRefGenome();
        refGenome.RefGenomeMap.put(CHR_1, refBases);
        ImpactClassifier classifier = new ImpactClassifier(refGenome);

        TranscriptData posTrans = createTransExons(
                GENE_ID_1, TRANS_ID_1, FusionCommon.POS_STRAND, new int[] { 20, 44, 68, 92 }, 11,
                20, 104, false, "");

        VariantData var = createDupVariant(refBases, 20, 3);
        VariantTransImpact impact = classifier.classifyVariant(var, posTrans);
        assertEquals("c.2_4dupTAG", impact.codingContext().Hgvs);

        var = createDupVariant(refBases, 31, 4);
        impact = classifier.classifyVariant(var, posTrans);
        assertEquals("c.12+1_12+4dupCGAG", impact.codingContext().Hgvs);

        var = createDupVariant(refBases, 33, 4);
        impact = classifier.classifyVariant(var, posTrans);
        assertEquals("c.12+3_12+6dupAGGG", impact.codingContext().Hgvs);

        //CGAGGGCCCAAA
        var = createDupVariant(refBases, 34, 4);
        impact = classifier.classifyVariant(var, posTrans);
        assertEquals("c.12+4_12+7dupGGGC", impact.codingContext().Hgvs); // should this be 12+4_13-5?

        var = createDupVariant(refBases, 36, 4);
        impact = classifier.classifyVariant(var, posTrans);
        assertEquals("c.12+6_12+9dupGCCC", impact.codingContext().Hgvs);

        var = createDupVariant(refBases, 37, 4);
        impact = classifier.classifyVariant(var, posTrans);
//        assertEquals("c.13-5_13-2dupCCAA", impact.codingContext().Hgvs);

        var = createDupVariant(refBases, 38, 4);
        impact = classifier.classifyVariant(var, posTrans);
        assertEquals("c.13-5_13-2dupCCAA", impact.codingContext().Hgvs);

        var = createDupVariant(refBases, 39, 4);
        impact = classifier.classifyVariant(var, posTrans);
        assertEquals("c.13-4_13-1dupCAAA", impact.codingContext().Hgvs);

        var = createDupVariant(refBases, 40, 4);
        impact = classifier.classifyVariant(var, posTrans);
        assertEquals("c.13-3_13-0dupAAAG", impact.codingContext().Hgvs);

        var = createDupVariant(refBases, 42, 2);
        impact = classifier.classifyVariant(var, posTrans);
//                assertEquals("c.13-6_11-3dupCCAA", impact.codingContext().Hgvs);
    }

    private VariantData createDupVariant(String refBases, int pos, int length)
    {
        String ref = refBases.substring(pos, pos + 1);
        String alt = refBases.substring(pos, pos + length + 1);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);
        String altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);
        return var;
    }
}
