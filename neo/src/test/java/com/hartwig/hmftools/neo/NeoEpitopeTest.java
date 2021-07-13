package com.hartwig.hmftools.neo;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.neo.NeoEpitopeType.STOP_LOST;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.gene.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.codon.AminoAcidConverter.reverseStrandBases;
import static com.hartwig.hmftools.neo.NeoEpitopeUtilsTest.CHR_1;
import static com.hartwig.hmftools.neo.NeoEpitopeUtilsTest.GENE_ID_1;
import static com.hartwig.hmftools.neo.NeoEpitopeUtilsTest.TRANS_ID_1;
import static com.hartwig.hmftools.neo.NeoEpitopeUtilsTest.generateRandomBases;

import static org.junit.Assert.assertEquals;

import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.MockRefGenome;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;

import org.junit.Assert;
import org.junit.Test;

public class NeoEpitopeTest
{
    @Test
    public void testSnvPointMutations()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateRandomBases(120);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        // tests: SNP, delete, insertion

        int[] exonStarts = {0, 20, 40, 60, 80, 100};
        Integer codingStart = new Integer(25);
        Integer codingEnd = new Integer(85);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        PointMutationData pmData = new PointMutationData(
                CHR_1, 42, chr1Bases.substring(42, 43), "A", GENE_ID_1, MISSENSE, 1, -1);

        NeoEpitope neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(42, neData.position(FS_UP));
        assertEquals(43, neData.position(FS_DOWN));
        assertEquals(POS_ORIENT, neData.orientation(FS_UP));
        assertEquals(NEG_ORIENT, neData.orientation(FS_DOWN));

        assertEquals(CODING, neData.CodingType[FS_UP]);
        assertEquals(CODING, neData.CodingType[FS_DOWN]);
        assertEquals(EXONIC, neData.RegionType[FS_UP]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWN]);
        assertEquals(PHASE_0, neData.Phases[FS_UP]);
        assertEquals(PHASE_1, neData.Phases[FS_DOWN]);
        assertEquals(3, neData.ExonRank[FS_UP]);
        assertEquals(3, neData.ExonRank[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        String upBases = chr1Bases.substring(25, 31);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        String novelBases = chr1Bases.substring(40, 42) + pmData.Alt;
        assertEquals(novelBases, neData.NovelCodonBases);

        String downBases = chr1Bases.substring(43, 51) + chr1Bases.substring(60, 61);
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);

        // SNP requiring all bases plus for phasing
        pmData = new PointMutationData(
                CHR_1, 46, chr1Bases.substring(46, 47), "A", GENE_ID_1, MISSENSE, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(46, neData.position(FS_UP));
        assertEquals(47, neData.position(FS_DOWN));
        assertEquals(POS_ORIENT, neData.orientation(FS_UP));
        assertEquals(NEG_ORIENT, neData.orientation(FS_DOWN));

        assertEquals(CODING, neData.CodingType[FS_UP]);
        assertEquals(CODING, neData.CodingType[FS_DOWN]);
        assertEquals(EXONIC, neData.RegionType[FS_UP]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWN]);
        assertEquals(PHASE_1, neData.Phases[FS_UP]);
        assertEquals(PHASE_2, neData.Phases[FS_DOWN]);
        assertEquals(3, neData.ExonRank[FS_UP]);
        assertEquals(3, neData.ExonRank[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(28, 31) + chr1Bases.substring(40, 46);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = pmData.Alt + chr1Bases.substring(47, 49);
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = chr1Bases.substring(49, 51) + chr1Bases.substring(60, 67);
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);

        // point mutation on last base of exon
        pmData = new PointMutationData(
                CHR_1, 30, chr1Bases.substring(30, 31), "A", GENE_ID_1, MISSENSE, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(30, neData.position(FS_UP));
        assertEquals(31, neData.position(FS_DOWN));
        assertEquals(POS_ORIENT, neData.orientation(FS_UP));
        assertEquals(NEG_ORIENT, neData.orientation(FS_DOWN));

        assertEquals(CODING, neData.CodingType[FS_UP]);
        assertEquals(CODING, neData.CodingType[FS_DOWN]);
        assertEquals(EXONIC, neData.RegionType[FS_UP]);
        assertEquals(INTRONIC, neData.RegionType[FS_DOWN]);
        assertEquals(PHASE_0, neData.Phases[FS_UP]);
        assertEquals(PHASE_0, neData.Phases[FS_DOWN]);
        assertEquals(2, neData.ExonRank[FS_UP]);
        assertEquals(3, neData.ExonRank[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(25, 28);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = chr1Bases.substring(28, 30) + pmData.Alt;
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = chr1Bases.substring(40, 49);
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);


        // test again on the reverse strand
        codingStart = 25;
        codingEnd = 65;
        TranscriptData transDataNegStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        pmData = new PointMutationData(
                CHR_1, 49, chr1Bases.substring(49, 50), "A", GENE_ID_1, MISSENSE, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataNegStrand, transDataNegStrand);

        assertEquals(49, neData.position(FS_UP));
        assertEquals(48, neData.position(FS_DOWN));
        assertEquals(NEG_ORIENT, neData.orientation(FS_UP));
        assertEquals(POS_ORIENT, neData.orientation(FS_DOWN));

        assertEquals(PHASE_2, neData.Phases[FS_UP]);
        assertEquals(PHASE_0, neData.Phases[FS_DOWN]);
        assertEquals(4, neData.ExonRank[FS_UP]);
        assertEquals(4, neData.ExonRank[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = reverseStrandBases(chr1Bases.substring(60, 66));
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = reverseStrandBases(chr1Bases.substring(48, 49) + pmData.Alt + chr1Bases.substring(50, 51));
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = reverseStrandBases(chr1Bases.substring(30, 31) + chr1Bases.substring(40, 48));
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);

        // point mutation on last base of exon
        pmData = new PointMutationData(
                CHR_1, 40, chr1Bases.substring(40, 41), "A", GENE_ID_1, MISSENSE, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataNegStrand, transDataNegStrand);

        assertEquals(40, neData.position(FS_UP));
        assertEquals(39, neData.position(FS_DOWN));
        assertEquals(EXONIC, neData.RegionType[FS_UP]);
        assertEquals(INTRONIC, neData.RegionType[FS_DOWN]);
        assertEquals(PHASE_2, neData.Phases[FS_UP]);
        assertEquals(PHASE_2, neData.Phases[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = reverseStrandBases(chr1Bases.substring(42, 51));
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = reverseStrandBases(chr1Bases.substring(30, 31) + pmData.Alt + chr1Bases.substring(41, 42));
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = reverseStrandBases(chr1Bases.substring(25, 30));
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);
    }

    @Test
    public void testMnvPointMutations()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateRandomBases(120);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        // tests: SNP, delete, insertion

        int[] exonStarts = {0, 20, 40, 60, 80, 100};
        Integer codingStart = new Integer(25);
        Integer codingEnd = new Integer(85);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        // mutation of 2 bases at first base of codon
        PointMutationData pmData = new PointMutationData(
                CHR_1, 43, chr1Bases.substring(43, 45), "GG", GENE_ID_1, MISSENSE, 1, -1);

        NeoEpitope neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(43, neData.position(FS_UP));
        assertEquals(45, neData.position(FS_DOWN));

        assertEquals(PHASE_1, neData.Phases[FS_UP]);
        assertEquals(PHASE_0, neData.Phases[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        String upBases = chr1Bases.substring(25, 31) + chr1Bases.substring(40, 43);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        String novelBases = pmData.Alt + chr1Bases.substring(45, 46);
        assertEquals(novelBases, neData.NovelCodonBases);

        String downBases = chr1Bases.substring(46, 51) + chr1Bases.substring(60, 64);
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);

        // 2-base MNV mutation on last base of exon, making 1 novel codon
        pmData = new PointMutationData(
                CHR_1, 30, chr1Bases.substring(30, 32),"AA", GENE_ID_1, MISSENSE, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(30, neData.position(FS_UP));
        assertEquals(32, neData.position(FS_DOWN));
        assertEquals(PHASE_0, neData.Phases[FS_UP]);
        assertEquals(PHASE_0, neData.Phases[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(25, 28);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        // FIXME:
        novelBases = chr1Bases.substring(28, 30) + pmData.Alt.substring(0, 1);
        // assertEquals(novelBases, neData.NovelCodonBases);

        downBases = chr1Bases.substring(40, 49);
        //assertEquals(downBases, neData.CodingBases[FS_DOWN]);


        // 5-base MNV mutation on second base of exon, making 3 novel codons
        pmData = new PointMutationData(
                CHR_1, 42, chr1Bases.substring(42, 47),
                "AAAAA", GENE_ID_1, MISSENSE, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(42, neData.position(FS_UP));
        assertEquals(47, neData.position(FS_DOWN));
        assertEquals(PHASE_0, neData.Phases[FS_UP]);
        assertEquals(PHASE_2, neData.Phases[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(25, 31);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = chr1Bases.substring(40, 42) + pmData.Alt + chr1Bases.substring(47, 49);
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = chr1Bases.substring(49, 51) + chr1Bases.substring(60, 67);
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);

        /*
        // test again on the reverse strand
        codingStart = 25;
        codingEnd = 65;
        TranscriptData transDataNegStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        pmData = new PointMutationData(
                CHR_1, 49, chr1Bases.substring(49, 50), "A", GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataNegStrand, transDataNegStrand);

        assertEquals(49, neData.position(FS_UP));
        assertEquals(48, neData.position(FS_DOWN));
        assertEquals(NEG_ORIENT, neData.orientation(FS_UP));
        assertEquals(POS_ORIENT, neData.orientation(FS_DOWN));

        assertEquals(CODING, neData.CodingType[FS_UP]);
        assertEquals(CODING, neData.CodingType[FS_DOWN]);
        assertEquals(EXONIC, neData.RegionType[FS_UP]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWN]);
        assertEquals(PHASE_2, neData.Phases[FS_UP]);
        assertEquals(PHASE_0, neData.Phases[FS_DOWN]);
        assertEquals(4, neData.ExonRank[FS_UP]);
        assertEquals(4, neData.ExonRank[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = reverseStrandBases(chr1Bases.substring(60, 66));
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = reverseStrandBases(chr1Bases.substring(48, 49) + pmData.Alt + chr1Bases.substring(50, 51));
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = reverseStrandBases(chr1Bases.substring(30, 31) + chr1Bases.substring(40, 48));
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);

        // point mutation on last base of exon
        pmData = new PointMutationData(
                CHR_1, 40, chr1Bases.substring(40, 41), "A", GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataNegStrand, transDataNegStrand);

        assertEquals(40, neData.position(FS_UP));
        assertEquals(39, neData.position(FS_DOWN));
        assertEquals(NEG_ORIENT, neData.orientation(FS_UP));
        assertEquals(POS_ORIENT, neData.orientation(FS_DOWN));

        assertEquals(CODING, neData.CodingType[FS_UP]);
        assertEquals(CODING, neData.CodingType[FS_DOWN]);
        assertEquals(EXONIC, neData.RegionType[FS_UP]);
        assertEquals(INTRONIC, neData.RegionType[FS_DOWN]);
        assertEquals(PHASE_2, neData.Phases[FS_UP]);
        assertEquals(PHASE_2, neData.Phases[FS_DOWN]);
        assertEquals(4, neData.ExonRank[FS_UP]);
        assertEquals(3, neData.ExonRank[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = reverseStrandBases(chr1Bases.substring(42, 51));
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = reverseStrandBases(chr1Bases.substring(30, 31) + pmData.Alt + chr1Bases.substring(41, 42));
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = reverseStrandBases(chr1Bases.substring(25, 30));
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);
        */
    }

    @Test
    public void testDeletionPointMutations()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateRandomBases(120);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        int[] exonStarts = { 0, 20, 40, 60, 80, 100 };
        Integer codingStart = new Integer(25);
        Integer codingEnd = new Integer(85);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        // delete - of 3 bases, in phase
        PointMutationData pmData = new PointMutationData(
                CHR_1, 41, chr1Bases.substring(41, 45), chr1Bases.substring(41, 42), GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        PmNeoEpitope neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(41, neData.position(FS_UP));
        assertEquals(45, neData.position(FS_DOWN));
        assertEquals(POS_ORIENT, neData.orientation(FS_UP));
        assertEquals(NEG_ORIENT, neData.orientation(FS_DOWN));

        assertEquals(CODING, neData.CodingType[FS_UP]);
        assertEquals(CODING, neData.CodingType[FS_DOWN]);
        assertEquals(EXONIC, neData.RegionType[FS_UP]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWN]);
        assertEquals(PHASE_2, neData.Phases[FS_UP]);
        assertEquals(PHASE_0, neData.Phases[FS_DOWN]);
        assertEquals(3, neData.ExonRank[FS_UP]);
        assertEquals(3, neData.ExonRank[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        String upBases = chr1Bases.substring(25, 31);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        String novelBases = chr1Bases.substring(40, 41) + pmData.Alt + chr1Bases.substring(45, 46);
        assertEquals(novelBases, neData.NovelCodonBases);

        String downBases = chr1Bases.substring(46, 51) + chr1Bases.substring(60, 64);
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);

        // delete - of 2 bases, out of phase
        pmData = new PointMutationData(
                CHR_1, 46, chr1Bases.substring(46, 49), chr1Bases.substring(46, 47),
                GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(46, neData.position(FS_UP));
        assertEquals(49, neData.position(FS_DOWN));
        assertEquals(PHASE_1, neData.Phases[FS_UP]);
        assertEquals(PHASE_1, neData.Phases[FS_DOWN]);
        Assert.assertFalse(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(28, 31) + chr1Bases.substring(40, 46);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = chr1Bases.substring(46, 47) + chr1Bases.substring(49, 51) + chr1Bases.substring(60, 71) + chr1Bases.substring(80, 91)
                + chr1Bases.substring(100, 111);
        assertEquals(novelBases, neData.NovelCodonBases);

        assertEquals("", neData.CodingBases[FS_DOWN]);

        // delete - of 2 bases, out of phase, starting mnid codon
        pmData = new PointMutationData(
                CHR_1, 47, chr1Bases.substring(47, 50), chr1Bases.substring(47, 48),
                GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(47, neData.position(FS_UP));
        assertEquals(50, neData.position(FS_DOWN));
        assertEquals(PHASE_2, neData.Phases[FS_UP]);
        assertEquals(PHASE_2, neData.Phases[FS_DOWN]);
        Assert.assertFalse(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(28, 31) + chr1Bases.substring(40, 46);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = chr1Bases.substring(46, 48) + chr1Bases.substring(50, 51) + chr1Bases.substring(60, 71) + chr1Bases.substring(80, 91)
                + chr1Bases.substring(100, 111);
        assertEquals(novelBases, neData.NovelCodonBases);

        assertEquals("", neData.CodingBases[FS_DOWN]);

        // delete - of 2 bases, out of phase, but deleting a new codon's first 2 bases
        pmData = new PointMutationData(
                CHR_1, 42, chr1Bases.substring(42, 45), chr1Bases.substring(42, 43),
                GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(42, neData.position(FS_UP));
        assertEquals(45, neData.position(FS_DOWN));
        assertEquals(PHASE_0, neData.Phases[FS_UP]);
        assertEquals(PHASE_0, neData.Phases[FS_DOWN]);
        Assert.assertFalse(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(25, 31) + chr1Bases.substring(40, 43);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = chr1Bases.substring(45, 51) + chr1Bases.substring(60, 71) + chr1Bases.substring(80, 91)
                + chr1Bases.substring(100, 111);
        assertEquals(novelBases, neData.NovelCodonBases);

        assertEquals("", neData.CodingBases[FS_DOWN]);

        // same again but with enough room to get required upstream bases
        pmData = new PointMutationData(
                CHR_1, 45, chr1Bases.substring(45, 48), chr1Bases.substring(45, 46),
                GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(45, neData.position(FS_UP));
        assertEquals(48, neData.position(FS_DOWN));
        assertEquals(PHASE_0, neData.Phases[FS_UP]);
        assertEquals(PHASE_0, neData.Phases[FS_DOWN]);
        Assert.assertFalse(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(28, 31) + chr1Bases.substring(40, 46);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = chr1Bases.substring(48, 51) + chr1Bases.substring(60, 71) + chr1Bases.substring(80, 91)
                + chr1Bases.substring(100, 111);
        assertEquals(novelBases, neData.NovelCodonBases);

        assertEquals("", neData.CodingBases[FS_DOWN]);


        // delete of 3 bases at last base of codon, no novel codon, in phase
        pmData = new PointMutationData(
                CHR_1, 42, chr1Bases.substring(42, 46), chr1Bases.substring(42, 43),
                GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(42, neData.position(FS_UP));
        assertEquals(46, neData.position(FS_DOWN));
        assertEquals(PHASE_0, neData.Phases[FS_UP]);
        assertEquals(PHASE_1, neData.Phases[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(25, 31) + chr1Bases.substring(40, 43);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        assertEquals("", neData.NovelCodonBases);

        downBases = chr1Bases.substring(46, 51) + chr1Bases.substring(60, 64);
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);

        neData.setAminoAcids(refGenome, 3);
        String upWildtypeBases = chr1Bases.substring(25, 31) + chr1Bases.substring(40, 51) + chr1Bases.substring(60, 61);
        String upWildAAs = NeoUtils.getAminoAcids(upWildtypeBases, true);
        assertEquals(upWildAAs, neData.UpstreamWildTypeAcids);

        // same again but allowing for enough coding bases upstream to satisfy the required count
        pmData = new PointMutationData(
                CHR_1, 45, chr1Bases.substring(45, 49), chr1Bases.substring(45, 46),
                GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(45, neData.position(FS_UP));
        assertEquals(49, neData.position(FS_DOWN));
        assertEquals(PHASE_0, neData.Phases[FS_UP]);
        assertEquals(PHASE_1, neData.Phases[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(28, 31) + chr1Bases.substring(40, 46);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        assertEquals("", neData.NovelCodonBases);

        downBases = chr1Bases.substring(49, 51) + chr1Bases.substring(60, 67);
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);


        // delete on negative strand
        codingStart = 25;
        codingEnd = 65;
        TranscriptData transDataNegStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        // 3 base DEL - deleting a codon
        pmData = new PointMutationData(
                CHR_1, 41, chr1Bases.substring(41, 45), chr1Bases.substring(41, 42), GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataNegStrand, transDataNegStrand);

        assertEquals(45, neData.position(FS_UP));
        assertEquals(41, neData.position(FS_DOWN));

        assertEquals(PHASE_0, neData.Phases[FS_UP]);
        assertEquals(PHASE_1, neData.Phases[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = reverseStrandBases(chr1Bases.substring(45, 51) + chr1Bases.substring(60, 63));
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        assertEquals("", neData.NovelCodonBases);

        downBases = reverseStrandBases(chr1Bases.substring(25, 31) + chr1Bases.substring(40, 42));
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);

        // 3 base DEL - crossing a codon boundary
        pmData = new PointMutationData(
                CHR_1, 42, chr1Bases.substring(42, 46), chr1Bases.substring(42, 43),
                GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataNegStrand, transDataNegStrand);

        assertEquals(46, neData.position(FS_UP));
        assertEquals(42, neData.position(FS_DOWN));

        assertEquals(PHASE_2, neData.Phases[FS_UP]);
        assertEquals(PHASE_0, neData.Phases[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = reverseStrandBases(chr1Bases.substring(48, 51) + chr1Bases.substring(60, 66));
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = reverseStrandBases(chr1Bases.substring(42, 43) + chr1Bases.substring(46, 48));
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = reverseStrandBases(chr1Bases.substring(25, 31) + chr1Bases.substring(40, 42));
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);

        neData.setAminoAcids(refGenome, 3);
        upWildtypeBases = chr1Bases.substring(30, 31) + chr1Bases.substring(40, 51) + chr1Bases.substring(60, 66);
        upWildAAs = NeoUtils.getAminoAcids(reverseStrandBases(upWildtypeBases), true);
        assertEquals(upWildAAs, neData.UpstreamWildTypeAcids);

        // again but a base up
        pmData = new PointMutationData(
                CHR_1, 43, chr1Bases.substring(43, 47), chr1Bases.substring(43, 44),
                GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataNegStrand, transDataNegStrand);

        assertEquals(47, neData.position(FS_UP));
        assertEquals(43, neData.position(FS_DOWN));

        assertEquals(PHASE_1, neData.Phases[FS_UP]);
        assertEquals(PHASE_2, neData.Phases[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = reverseStrandBases(chr1Bases.substring(48, 51) + chr1Bases.substring(60, 66));
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = reverseStrandBases(chr1Bases.substring(42, 44) + chr1Bases.substring(47, 48));
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = reverseStrandBases(chr1Bases.substring(25, 31) + chr1Bases.substring(40, 42));
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);
    }

    @Test
    public void testInsertPointMutations()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateRandomBases(120);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        int[] exonStarts = { 0, 20, 40, 60, 80, 100 };
        Integer codingStart = new Integer(25);
        Integer codingEnd = new Integer(85);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        // insert of 2 bases at first base of codon - frameshift
        PointMutationData pmData = new PointMutationData(
                CHR_1, 28, chr1Bases.substring(28, 29), "AAA", GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        NeoEpitope neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(28, neData.position(FS_UP));
        assertEquals(29, neData.position(FS_DOWN));
        assertEquals(POS_ORIENT, neData.orientation(FS_UP));
        assertEquals(NEG_ORIENT, neData.orientation(FS_DOWN));

        assertEquals(CODING, neData.CodingType[FS_UP]);
        assertEquals(CODING, neData.CodingType[FS_DOWN]);
        assertEquals(EXONIC, neData.RegionType[FS_UP]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWN]);
        assertEquals(PHASE_1, neData.Phases[FS_UP]);
        assertEquals(PHASE_2, neData.Phases[FS_DOWN]); // whatever it would usually be
        assertEquals(2, neData.ExonRank[FS_UP]);
        assertEquals(2, neData.ExonRank[FS_DOWN]);
        assertFalse(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        String upBases = chr1Bases.substring(25, 28);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        String downBases = pmData.Alt + chr1Bases.substring(29, 31) + chr1Bases.substring(40, 51) + chr1Bases.substring(60, 71) + chr1Bases.substring(80, 91)
                + chr1Bases.substring(100, 111);

        assertEquals("", neData.CodingBases[FS_DOWN]);
        assertEquals(downBases, neData.NovelCodonBases); // due to lack of phasing

        // insert of 3 bases at last base of codon - making 1 novel codon, starting after first base and in phase
        pmData = new PointMutationData(
                CHR_1, 42, chr1Bases.substring(42, 43), chr1Bases.substring(42, 43) + "AAA",
                GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(42, neData.position(FS_UP));
        assertEquals(43, neData.position(FS_DOWN));
        assertEquals(PHASE_0, neData.Phases[FS_UP]); // pushed out 3 by the inserted bases, unch
        assertEquals(PHASE_1, neData.Phases[FS_DOWN]);
        assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(25, 31) + chr1Bases.substring(40, 43);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        String novelBases = pmData.Alt.substring(1);
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = chr1Bases.substring(43, 51) + chr1Bases.substring(60, 61);

        assertEquals(downBases, neData.CodingBases[FS_DOWN]);

        // insert of 3 bases at 2nd base of codon - 2 novel codons, starting after first base and in phase
        pmData = new PointMutationData(
                CHR_1, 41, chr1Bases.substring(41, 42), chr1Bases.substring(41, 42) + "AAA",
                GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(41, neData.position(FS_UP));
        assertEquals(42, neData.position(FS_DOWN));
        assertEquals(PHASE_2, neData.Phases[FS_UP]); // pushed out 6 by the inserted bases, unch
        assertEquals(PHASE_0, neData.Phases[FS_DOWN]);
        assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(25, 31);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = chr1Bases.substring(40, 41) + pmData.Alt + chr1Bases.substring(42, 43);
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = chr1Bases.substring(43, 51) + chr1Bases.substring(60, 61);

        assertEquals(downBases, neData.CodingBases[FS_DOWN]);

        // negative strand

        // insert of 3 bases, in phase
        codingStart = 25;
        codingEnd = 65;
        TranscriptData transDataNegStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        pmData = new PointMutationData(
                CHR_1, 43, chr1Bases.substring(43, 44), "AAAA", GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataNegStrand, transDataNegStrand);

        assertEquals(44, neData.position(FS_UP));
        assertEquals(43, neData.position(FS_DOWN));
        assertEquals(PHASE_1, neData.Phases[FS_UP]); // pushed out 2 by the inserted bases
        assertEquals(PHASE_2, neData.Phases[FS_DOWN]); // whatever it would usually be
        assertEquals(4, neData.ExonRank[FS_UP]);
        assertEquals(4, neData.ExonRank[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = reverseStrandBases(chr1Bases.substring(45, 51) + chr1Bases.substring(60, 63));
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = reverseStrandBases(chr1Bases.substring(42, 43) + pmData.Alt + chr1Bases.substring(44, 45));
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = reverseStrandBases(chr1Bases.substring(25, 31) + chr1Bases.substring(40, 42));
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);

        // insert of 3 bases, in phase, with inserted bases forming the new codon
        pmData = new PointMutationData(
                CHR_1, 44, chr1Bases.substring(44, 45), chr1Bases.substring(44, 45) + "AAA",
                GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataNegStrand, transDataNegStrand);

        assertEquals(45, neData.position(FS_UP));
        assertEquals(44, neData.position(FS_DOWN));
        assertEquals(PHASE_0, neData.Phases[FS_UP]); // pushed out 2 by the inserted bases
        assertEquals(PHASE_1, neData.Phases[FS_DOWN]); // whatever it would usually be
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = reverseStrandBases(chr1Bases.substring(45, 51) + chr1Bases.substring(60, 63));
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = reverseStrandBases(chr1Bases.substring(42, 44) + pmData.Alt);
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = reverseStrandBases(chr1Bases.substring(25, 31) + chr1Bases.substring(40, 42));
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);
    }

    @Test
    public void testStopLostPointMutations()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateRandomBases(64) + "TGA" + generateRandomBases(60);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        int[] exonStarts = { 0, 20, 40, 60, 80, 100 };
        Integer codingStart = new Integer(25);
        Integer codingEnd = new Integer(66);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        // stop-lost
        PointMutationData pmData = new PointMutationData(
                CHR_1, 64, chr1Bases.substring(64, 65), "A", GENE_ID_1, MISSENSE, 1, -1);

        PmNeoEpitope neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(CODING, neData.CodingType[FS_UP]);
        assertEquals(CODING, neData.CodingType[FS_DOWN]);
        assertEquals(EXONIC, neData.RegionType[FS_UP]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWN]);
        assertEquals(PHASE_1, neData.Phases[FS_UP]);
        assertEquals(PHASE_2, neData.Phases[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        assertEquals(NeoEpitopeType.MISSENSE, neData.variantType());

        neData.setCodingBases(refGenome, 3);
        neData.setAminoAcids(refGenome, 3);
        assertEquals(STOP_LOST, neData.variantType());

        String upBases = chr1Bases.substring(46, 51) + chr1Bases.substring(60, 64);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        String novelBases = pmData.Alt + chr1Bases.substring(65, 71) + chr1Bases.substring(80, 91) + chr1Bases.substring(100, 111);
        assertEquals(novelBases, neData.NovelCodonBases);

        assertEquals("", neData.CodingBases[FS_DOWN]);
    }

    @Test
    public void testSvFusions()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateRandomBases(120);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        // tests: intronic vs exonic, with and without insertions

        int[] exonStarts = { 0, 20, 40, 60, 80, 100 };
        Integer codingStart = new Integer(25);
        Integer codingEnd = new Integer(85);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        final String[] validTrans = new String[] {transDataPosStrand.TransName, transDataPosStrand.TransName};

        // intronic
        NeoEpitopeFusion svData = new NeoEpitopeFusion(
                GENE_ID_1, GENE_ID_1, CHR_1, 35, POS_ORIENT, 1, GENE_ID_1, GENE_ID_1, CHR_1, 55, NEG_ORIENT,
                1, 1, "", 0, validTrans);

        NeoEpitope neData = new SvNeoEpitope(svData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(35, neData.position(FS_UP));
        assertEquals(55, neData.position(FS_DOWN));
        assertEquals(POS_ORIENT, neData.orientation(FS_UP));
        assertEquals(NEG_ORIENT, neData.orientation(FS_DOWN));

        assertEquals(CODING, neData.CodingType[FS_UP]);
        assertEquals(CODING, neData.CodingType[FS_DOWN]);
        assertEquals(INTRONIC, neData.RegionType[FS_UP]);
        assertEquals(INTRONIC, neData.RegionType[FS_DOWN]);
        assertEquals(PHASE_0, neData.Phases[FS_UP]);
        assertEquals(PHASE_2, neData.Phases[FS_DOWN]);
        assertEquals(2, neData.ExonRank[FS_UP]);
        assertEquals(4, neData.ExonRank[FS_DOWN]);
        assertFalse(neData.phaseMatched());

        neData.setCodingBases(refGenome, 2);

        String upBases = chr1Bases.substring(25, 31);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        String novelBases = chr1Bases.substring(60, 71) + chr1Bases.substring(80, 91) + chr1Bases.substring(100, 111);

        assertEquals(novelBases, neData.NovelCodonBases);

        assertEquals("", neData.CodingBases[FS_DOWN]);

        neData.setAminoAcids(refGenome, 3);
        String upWildtypeBases = chr1Bases.substring(25, 31) + chr1Bases.substring(40, 50);
        String upWildAAs = NeoUtils.getAminoAcids(upWildtypeBases, true);
        assertEquals(upWildAAs, neData.UpstreamWildTypeAcids);

        // intronic to exonic - skips to next exon
        svData = new NeoEpitopeFusion(
                GENE_ID_1, GENE_ID_1, CHR_1, 35, POS_ORIENT, 1, GENE_ID_1, GENE_ID_1, CHR_1, 45, NEG_ORIENT,
                1, 1, "", 0, validTrans);

        neData = new SvNeoEpitope(svData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(CODING, neData.CodingType[FS_UP]);
        assertEquals(CODING, neData.CodingType[FS_DOWN]);
        assertEquals(INTRONIC, neData.RegionType[FS_UP]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWN]);
        assertEquals(PHASE_0, neData.Phases[FS_UP]);
        assertEquals(PHASE_2, neData.Phases[FS_DOWN]);
        assertEquals(2, neData.ExonRank[FS_UP]);
        assertEquals(4, neData.ExonRank[FS_DOWN]);
        assertFalse(neData.phaseMatched());

        neData.setCodingBases(refGenome, 2);

        upBases = chr1Bases.substring(25, 31);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = chr1Bases.substring(60, 71) + chr1Bases.substring(80, 91) + chr1Bases.substring(100, 111);
        assertEquals(novelBases, neData.NovelCodonBases);

        assertEquals("", neData.CodingBases[FS_DOWN]);

        // exonic to exonic, in phase
        svData = new NeoEpitopeFusion(
                GENE_ID_1, GENE_ID_1, CHR_1, 44, POS_ORIENT, 1, GENE_ID_1, GENE_ID_1, CHR_1, 63, NEG_ORIENT,
                1, 1, "", 0, validTrans);

        neData = new SvNeoEpitope(svData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(CODING, neData.CodingType[FS_UP]);
        assertEquals(CODING, neData.CodingType[FS_DOWN]);
        assertEquals(EXONIC, neData.RegionType[FS_UP]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWN]);
        assertEquals(PHASE_2, neData.Phases[FS_UP]);
        assertEquals(PHASE_0, neData.Phases[FS_DOWN]);
        assertEquals(3, neData.ExonRank[FS_UP]);
        assertEquals(4, neData.ExonRank[FS_DOWN]);
        Assert.assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 2);

        upBases = chr1Bases.substring(28, 31) + chr1Bases.substring(40, 43);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = chr1Bases.substring(43, 45) + chr1Bases.substring(63, 64);
        assertEquals(novelBases, neData.NovelCodonBases);

        String downBases = chr1Bases.substring(64, 70);
        assertEquals(downBases, neData.CodingBases[FS_DOWN]);

        // exonic to exonic with an insert sequence, out of phase
        String insSequence = "AA";

        svData = new NeoEpitopeFusion(
                GENE_ID_1, GENE_ID_1, CHR_1, 44, POS_ORIENT, 1, GENE_ID_1, GENE_ID_1, CHR_1, 63, NEG_ORIENT,
                1, 1, insSequence, 0, validTrans);

        neData = new SvNeoEpitope(svData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(CODING, neData.CodingType[FS_UP]);
        assertEquals(CODING, neData.CodingType[FS_DOWN]);
        assertEquals(EXONIC, neData.RegionType[FS_UP]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWN]);
        assertEquals(PHASE_1, neData.Phases[FS_UP]);
        assertEquals(PHASE_0, neData.Phases[FS_DOWN]);
        assertEquals(3, neData.ExonRank[FS_UP]);
        assertEquals(4, neData.ExonRank[FS_DOWN]);
        assertFalse(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(25, 31) + chr1Bases.substring(40, 43);
        assertEquals(upBases, neData.CodingBases[FS_UP]);

        novelBases = chr1Bases.substring(43, 45) + insSequence + chr1Bases.substring(63, 71)
                + chr1Bases.substring(80, 91) + chr1Bases.substring(100, 111);
        assertEquals(novelBases, neData.NovelCodonBases);

        assertEquals("", neData.CodingBases[FS_DOWN]);

        neData.setAminoAcids(refGenome, 3);
        upWildtypeBases = chr1Bases.substring(25, 31) + chr1Bases.substring(40, 51) + chr1Bases.substring(60, 61);
        upWildAAs = NeoUtils.getAminoAcids(upWildtypeBases, true);
        assertEquals(upWildAAs, neData.UpstreamWildTypeAcids);
    }

}
