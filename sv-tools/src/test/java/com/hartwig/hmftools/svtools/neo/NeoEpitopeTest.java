package com.hartwig.hmftools.svtools.neo;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.reverseStrandBases;
import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.swapDnaToRna;
import static com.hartwig.hmftools.svtools.neo.AminoAcidConverter.swapRnaToDna;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.getCodingBases;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertFalse;

import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.MockRefGenome;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;

import org.junit.Test;

public class NeoEpitopeTest
{
    private static final String CHR_1 = "";
    private static final String GENE_ID_1 = "ENSG001";
    private static final int TRANS_ID_1 = 1;

    @Test
    public void testDnaRnaRoutines()
    {
        String dnaBases = "AGCT";
        String rnaBases = swapDnaToRna(dnaBases);
        assertTrue(rnaBases.equals("AGCU"));
        assertTrue(dnaBases.equals(swapRnaToDna(rnaBases)));

        dnaBases = "AGCTTCGACT";
        String reverseStrandDna = reverseStrandBases(dnaBases);
        assertTrue(reverseStrandDna.equals("AGTCGAAGCT"));
    }

    public static String generateRandomBases(int length)
    {
        char[] str = new char[length];
        String bases = "ACGT";

        int baseIndex = 0;
        for(int i = 0; i < length; ++i)
        {
            str[i] = bases.charAt(baseIndex);

            if(baseIndex == 3)
                baseIndex = 0;
            else
                ++baseIndex;
        }

        return String.valueOf(str);
    }

    @Test
    public void testSetCodingBases()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateRandomBases(100);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        // tests: repeated for each strand
        // upstream: intronic
        // upstream: intronic, coding bases insufficient for required bases
        // upstream: exonic
        // upstream: exonic, coding bases insufficient for required bases
        // upstream: exonic, coding finishes in same exon

        int[] exonStarts = {0, 20, 40, 60, 80, 100};
        Integer codingStart = new Integer(25);
        Integer codingEnd = new Integer(85);

        TranscriptData transDataUp = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        // upstream: intronic
        int nePosition = 32; // intronic
        byte neOrientation = NEG_ORIENT;

        String codingBases = getCodingBases(refGenome, transDataUp, CHR_1, nePosition, neOrientation, 9, true);

        String actCodingBases = chr1Bases.substring(40, 49);
        assertEquals(actCodingBases, codingBases);

        // spanning 3 exons
        codingBases = getCodingBases(refGenome, transDataUp, CHR_1, nePosition, neOrientation, 25, true);

        actCodingBases = chr1Bases.substring(40, 51) + chr1Bases.substring(60, 71) + chr1Bases.substring(80, 83);
        assertEquals(actCodingBases, codingBases);

        // upstream: intronic, coding bases insufficient for required bases
        nePosition = 76;
        codingBases = getCodingBases(refGenome, transDataUp, CHR_1, nePosition, neOrientation, 10, true);

        actCodingBases = chr1Bases.substring(80, 86);
        assertEquals(actCodingBases, codingBases);

        // upstream: exonic
        nePosition = 44;
        codingBases = getCodingBases(refGenome, transDataUp, CHR_1, nePosition, neOrientation, 10, true);

        actCodingBases = chr1Bases.substring(44, 51) + chr1Bases.substring(60, 63);
        assertEquals(actCodingBases, codingBases);

        // upstream: exonic, coding bases insufficient for required bases
        nePosition = 68;
        codingBases = getCodingBases(refGenome, transDataUp, CHR_1, nePosition, neOrientation, 10, true);

        actCodingBases = chr1Bases.substring(68, 71) + chr1Bases.substring(80, 86);
        assertEquals(actCodingBases, codingBases);

        // upstream: exonic, coding finishes in same exon
        nePosition = 81;
        codingBases = getCodingBases(refGenome, transDataUp, CHR_1, nePosition, neOrientation, 10, true);

        actCodingBases = chr1Bases.substring(81, 86);
        assertEquals(actCodingBases, codingBases);

        // downstream, when not allowed to start in the same exon
        nePosition = 44;
        codingBases = getCodingBases(refGenome, transDataUp, CHR_1, nePosition, neOrientation, 15, false);

        actCodingBases = chr1Bases.substring(60, 71) + chr1Bases.substring(80, 84);
        assertEquals(actCodingBases, codingBases);


        // test again with reverse orientation
        neOrientation = POS_ORIENT;
        nePosition = 55; // intronic

        codingBases = getCodingBases(refGenome, transDataUp, CHR_1, nePosition, neOrientation, 9, true);

        actCodingBases = chr1Bases.substring(42, 51);
        assertEquals(actCodingBases, codingBases);

        // spanning 3 exons
        nePosition = 75;
        codingBases = getCodingBases(refGenome, transDataUp, CHR_1, nePosition, neOrientation, 25, true);

        actCodingBases = chr1Bases.substring(28, 31) + chr1Bases.substring(40, 51) + chr1Bases.substring(60, 71);
        assertEquals(actCodingBases, codingBases);

        // downstream: intronic, coding bases insufficient for required bases
        nePosition = 35;
        codingBases = getCodingBases(refGenome, transDataUp, CHR_1, nePosition, neOrientation, 10, true);

        actCodingBases = chr1Bases.substring(25, 31);
        assertEquals(actCodingBases, codingBases);

        // downstream: exonic
        nePosition = 44;
        codingBases = getCodingBases(refGenome, transDataUp, CHR_1, nePosition, neOrientation, 10, true);

        actCodingBases = chr1Bases.substring(26, 31) + chr1Bases.substring(40, 45);
        assertEquals(actCodingBases, codingBases);

        // downstream: exonic, coding bases insufficient for required bases
        nePosition = 40;
        codingBases = getCodingBases(refGenome, transDataUp, CHR_1, nePosition, neOrientation, 10, true);

        actCodingBases = chr1Bases.substring(25, 31) + chr1Bases.substring(40, 41);
        assertEquals(actCodingBases, codingBases);

        // downstream: exonic, coding finishes in same exon
        nePosition = 29;
        codingBases = getCodingBases(refGenome, transDataUp, CHR_1, nePosition, neOrientation, 10, true);

        actCodingBases = chr1Bases.substring(25, 30);
        assertEquals(actCodingBases, codingBases);

        // downstream, when not allowed to start in the same exon
        nePosition = 64;
        codingBases = getCodingBases(refGenome, transDataUp, CHR_1, nePosition, neOrientation, 15, false);

        actCodingBases = chr1Bases.substring(27, 31) + chr1Bases.substring(40, 51);
        assertEquals(actCodingBases, codingBases);
    }

    @Test
    public void testPointMutations()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateRandomBases(100);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        // tests: SNP, delete, insertion

        int[] exonStarts = {0, 20, 40, 60, 80, 100};
        Integer codingStart = new Integer(25);
        Integer codingEnd = new Integer(85);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        PointMutationData pmData = new PointMutationData(
                CHR_1, 42, chr1Bases.substring(42, 43), "A", GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        NeoEpitope neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(42, neData.position(FS_UPSTREAM));
        assertEquals(43, neData.position(FS_DOWNSTREAM));
        assertEquals(POS_ORIENT, neData.orientation(FS_UPSTREAM));
        assertEquals(NEG_ORIENT, neData.orientation(FS_DOWNSTREAM));

        assertEquals(CODING, neData.CodingType[FS_UPSTREAM]);
        assertEquals(CODING, neData.CodingType[FS_DOWNSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_UPSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWNSTREAM]);
        assertEquals(2, neData.Phases[FS_UPSTREAM]);
        assertEquals(0, neData.Phases[FS_DOWNSTREAM]);
        assertEquals(3, neData.ExonRank[FS_UPSTREAM]);
        assertEquals(3, neData.ExonRank[FS_DOWNSTREAM]);
        assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        String upBases = chr1Bases.substring(25, 31) + chr1Bases.substring(40, 42) + pmData.Alt;
        assertEquals(upBases, neData.CodingBases[FS_UPSTREAM]);

        String downBases = chr1Bases.substring(43, 51) + chr1Bases.substring(60, 61);
        assertEquals(downBases, neData.CodingBases[FS_DOWNSTREAM]);

        // SNP requiring all bases plus for phasing
        pmData = new PointMutationData(
                CHR_1, 46, chr1Bases.substring(46, 47), "A", GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(46, neData.position(FS_UPSTREAM));
        assertEquals(47, neData.position(FS_DOWNSTREAM));
        assertEquals(POS_ORIENT, neData.orientation(FS_UPSTREAM));
        assertEquals(NEG_ORIENT, neData.orientation(FS_DOWNSTREAM));

        assertEquals(CODING, neData.CodingType[FS_UPSTREAM]);
        assertEquals(CODING, neData.CodingType[FS_DOWNSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_UPSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWNSTREAM]);
        assertEquals(0, neData.Phases[FS_UPSTREAM]);
        assertEquals(1, neData.Phases[FS_DOWNSTREAM]);
        assertEquals(3, neData.ExonRank[FS_UPSTREAM]);
        assertEquals(3, neData.ExonRank[FS_DOWNSTREAM]);
        assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(28, 31) + chr1Bases.substring(40, 46);
        assertEquals(upBases, neData.CodingBases[FS_UPSTREAM]);

        String novelBases = pmData.Alt + chr1Bases.substring(47, 49);
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = chr1Bases.substring(49, 51) + chr1Bases.substring(60, 67);
        assertEquals(downBases, neData.CodingBases[FS_DOWNSTREAM]);

        // point mutation on last base of exon
        pmData = new PointMutationData(
                CHR_1, 30, chr1Bases.substring(30, 31), "A", GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(30, neData.position(FS_UPSTREAM));
        assertEquals(31, neData.position(FS_DOWNSTREAM));
        assertEquals(POS_ORIENT, neData.orientation(FS_UPSTREAM));
        assertEquals(NEG_ORIENT, neData.orientation(FS_DOWNSTREAM));

        assertEquals(CODING, neData.CodingType[FS_UPSTREAM]);
        assertEquals(CODING, neData.CodingType[FS_DOWNSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_UPSTREAM]);
        assertEquals(INTRONIC, neData.RegionType[FS_DOWNSTREAM]);
        assertEquals(2, neData.Phases[FS_UPSTREAM]);
        assertEquals(2, neData.Phases[FS_DOWNSTREAM]);
        assertEquals(2, neData.ExonRank[FS_UPSTREAM]);
        assertEquals(3, neData.ExonRank[FS_DOWNSTREAM]);
        assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(25, 30) + pmData.Alt;
        assertEquals(upBases, neData.CodingBases[FS_UPSTREAM]);

        downBases = chr1Bases.substring(40, 49);
        assertEquals(downBases, neData.CodingBases[FS_DOWNSTREAM]);

        // insert
        pmData = new PointMutationData(
                CHR_1, 28, chr1Bases.substring(28, 29), "AAA", GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(28, neData.position(FS_UPSTREAM));
        assertEquals(29, neData.position(FS_DOWNSTREAM));
        assertEquals(POS_ORIENT, neData.orientation(FS_UPSTREAM));
        assertEquals(NEG_ORIENT, neData.orientation(FS_DOWNSTREAM));

        assertEquals(CODING, neData.CodingType[FS_UPSTREAM]);
        assertEquals(CODING, neData.CodingType[FS_DOWNSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_UPSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWNSTREAM]);
        assertEquals(2, neData.Phases[FS_UPSTREAM]); // pushed out 2 by the inserted bases
        assertEquals(1, neData.Phases[FS_DOWNSTREAM]); // whatever it would usually be
        assertEquals(2, neData.ExonRank[FS_UPSTREAM]);
        assertEquals(2, neData.ExonRank[FS_DOWNSTREAM]);
        assertFalse(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(25, 28) + pmData.Alt;
        assertEquals(upBases, neData.CodingBases[FS_UPSTREAM]);

        downBases = chr1Bases.substring(29, 31) + chr1Bases.substring(40, 47);
        assertEquals("", neData.CodingBases[FS_DOWNSTREAM]);
        assertEquals(downBases, neData.NovelCodonBases); // due to lack of phasing

        // delete
        pmData = new PointMutationData(
                CHR_1, 41, chr1Bases.substring(41, 45), chr1Bases.substring(41, 42), GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(41, neData.position(FS_UPSTREAM));
        assertEquals(45, neData.position(FS_DOWNSTREAM));
        assertEquals(POS_ORIENT, neData.orientation(FS_UPSTREAM));
        assertEquals(NEG_ORIENT, neData.orientation(FS_DOWNSTREAM));

        assertEquals(CODING, neData.CodingType[FS_UPSTREAM]);
        assertEquals(CODING, neData.CodingType[FS_DOWNSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_UPSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWNSTREAM]);
        assertEquals(1, neData.Phases[FS_UPSTREAM]);
        assertEquals(2, neData.Phases[FS_DOWNSTREAM]);
        assertEquals(3, neData.ExonRank[FS_UPSTREAM]);
        assertEquals(3, neData.ExonRank[FS_DOWNSTREAM]);
        assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(25, 31);
        assertEquals(upBases, neData.CodingBases[FS_UPSTREAM]);

        novelBases = chr1Bases.substring(40, 41) + pmData.Alt + chr1Bases.substring(45, 46);
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = chr1Bases.substring(46, 51) + chr1Bases.substring(60, 64);
        assertEquals(downBases, neData.CodingBases[FS_DOWNSTREAM]);

        // test again on the reverse strand
        codingStart = 25;
        codingEnd = 65;
        TranscriptData transDataNegStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        pmData = new PointMutationData(
                CHR_1, 49, chr1Bases.substring(49, 50), "A", GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataNegStrand, transDataNegStrand);

        assertEquals(49, neData.position(FS_UPSTREAM));
        assertEquals(48, neData.position(FS_DOWNSTREAM));
        assertEquals(NEG_ORIENT, neData.orientation(FS_UPSTREAM));
        assertEquals(POS_ORIENT, neData.orientation(FS_DOWNSTREAM));

        assertEquals(CODING, neData.CodingType[FS_UPSTREAM]);
        assertEquals(CODING, neData.CodingType[FS_DOWNSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_UPSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWNSTREAM]);
        assertEquals(1, neData.Phases[FS_UPSTREAM]);
        assertEquals(2, neData.Phases[FS_DOWNSTREAM]);
        assertEquals(4, neData.ExonRank[FS_UPSTREAM]);
        assertEquals(4, neData.ExonRank[FS_DOWNSTREAM]);
        assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = reverseStrandBases(chr1Bases.substring(60, 66));
        assertEquals(upBases, neData.CodingBases[FS_UPSTREAM]);

        novelBases = reverseStrandBases(chr1Bases.substring(48, 49) + pmData.Alt + chr1Bases.substring(50, 51));
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = reverseStrandBases(chr1Bases.substring(30, 31) + chr1Bases.substring(40, 48));
        assertEquals(downBases, neData.CodingBases[FS_DOWNSTREAM]);

        // point mutation on last base of exon
        pmData = new PointMutationData(
                CHR_1, 40, chr1Bases.substring(40, 41), "A", GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataNegStrand, transDataNegStrand);

        assertEquals(40, neData.position(FS_UPSTREAM));
        assertEquals(39, neData.position(FS_DOWNSTREAM));
        assertEquals(NEG_ORIENT, neData.orientation(FS_UPSTREAM));
        assertEquals(POS_ORIENT, neData.orientation(FS_DOWNSTREAM));

        assertEquals(CODING, neData.CodingType[FS_UPSTREAM]);
        assertEquals(CODING, neData.CodingType[FS_DOWNSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_UPSTREAM]);
        assertEquals(INTRONIC, neData.RegionType[FS_DOWNSTREAM]);
        assertEquals(1, neData.Phases[FS_UPSTREAM]);
        assertEquals(1, neData.Phases[FS_DOWNSTREAM]);
        assertEquals(4, neData.ExonRank[FS_UPSTREAM]);
        assertEquals(3, neData.ExonRank[FS_DOWNSTREAM]);
        assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = reverseStrandBases(chr1Bases.substring(42, 51));
        assertEquals(upBases, neData.CodingBases[FS_UPSTREAM]);

        novelBases = reverseStrandBases(chr1Bases.substring(30, 31) + pmData.Alt + chr1Bases.substring(41, 42));
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = reverseStrandBases(chr1Bases.substring(25, 30));
        assertEquals(downBases, neData.CodingBases[FS_DOWNSTREAM]);

        // insert
        pmData = new PointMutationData(
                CHR_1, 43, chr1Bases.substring(43, 44), "AAAA", GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataNegStrand, transDataNegStrand);

        assertEquals(43, neData.position(FS_UPSTREAM));
        assertEquals(42, neData.position(FS_DOWNSTREAM));
        assertEquals(NEG_ORIENT, neData.orientation(FS_UPSTREAM));
        assertEquals(POS_ORIENT, neData.orientation(FS_DOWNSTREAM));

        assertEquals(CODING, neData.CodingType[FS_UPSTREAM]);
        assertEquals(CODING, neData.CodingType[FS_DOWNSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_UPSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWNSTREAM]);
        assertEquals(1, neData.Phases[FS_UPSTREAM]); // pushed out 2 by the inserted bases
        assertEquals(2, neData.Phases[FS_DOWNSTREAM]); // whatever it would usually be
        assertEquals(4, neData.ExonRank[FS_UPSTREAM]);
        assertEquals(4, neData.ExonRank[FS_DOWNSTREAM]);
        assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = reverseStrandBases(pmData.Alt.substring(2) + chr1Bases.substring(44, 51));
        assertEquals(upBases, neData.CodingBases[FS_UPSTREAM]);

        novelBases = reverseStrandBases(chr1Bases.substring(42, 43) + pmData.Alt.substring(0, 2));
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = reverseStrandBases(chr1Bases.substring(25, 31) + chr1Bases.substring(40, 42));
        assertEquals(downBases, neData.CodingBases[FS_DOWNSTREAM]);

        // delete
        pmData = new PointMutationData(
                CHR_1, 41, chr1Bases.substring(41, 45), chr1Bases.substring(41, 42), GENE_ID_1, NONSENSE_OR_FRAMESHIFT, 1, -1);

        neData = new PmNeoEpitope(pmData);

        neData.setTranscriptData(transDataNegStrand, transDataNegStrand);

        assertEquals(41, neData.position(FS_UPSTREAM));
        assertEquals(40, neData.position(FS_DOWNSTREAM));
        assertEquals(NEG_ORIENT, neData.orientation(FS_UPSTREAM));
        assertEquals(POS_ORIENT, neData.orientation(FS_DOWNSTREAM));

        assertEquals(CODING, neData.CodingType[FS_UPSTREAM]);
        assertEquals(CODING, neData.CodingType[FS_DOWNSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_UPSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWNSTREAM]);
        assertEquals(0, neData.Phases[FS_UPSTREAM]);
        assertEquals(1, neData.Phases[FS_DOWNSTREAM]);
        assertEquals(4, neData.ExonRank[FS_UPSTREAM]);
        assertEquals(4, neData.ExonRank[FS_DOWNSTREAM]);
        assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = reverseStrandBases(chr1Bases.substring(45, 51) + chr1Bases.substring(60, 63));
        assertEquals(upBases, neData.CodingBases[FS_UPSTREAM]);

        novelBases = reverseStrandBases(chr1Bases.substring(30, 31) + chr1Bases.substring(40, 42));
        assertEquals(novelBases, neData.NovelCodonBases);

        downBases = reverseStrandBases(chr1Bases.substring(25, 30));
        assertEquals(downBases, neData.CodingBases[FS_DOWNSTREAM]);
    }

    @Test
    public void testSvFusions()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateRandomBases(100);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        // tests: intronic vs exonic, with and without insertions

        int[] exonStarts = { 0, 20, 40, 60, 80, 100 };
        Integer codingStart = new Integer(25);
        Integer codingEnd = new Integer(85);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        // intronic
        NeoEpitopeFusion svData = new NeoEpitopeFusion(
                GENE_ID_1, GENE_ID_1, CHR_1, 35, POS_ORIENT, 1, GENE_ID_1, GENE_ID_1, CHR_1, 55, NEG_ORIENT,
                1, 1, "");

        NeoEpitope neData = new SvNeoEpitope(svData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(35, neData.position(FS_UPSTREAM));
        assertEquals(55, neData.position(FS_DOWNSTREAM));
        assertEquals(POS_ORIENT, neData.orientation(FS_UPSTREAM));
        assertEquals(NEG_ORIENT, neData.orientation(FS_DOWNSTREAM));

        assertEquals(CODING, neData.CodingType[FS_UPSTREAM]);
        assertEquals(CODING, neData.CodingType[FS_DOWNSTREAM]);
        assertEquals(INTRONIC, neData.RegionType[FS_UPSTREAM]);
        assertEquals(INTRONIC, neData.RegionType[FS_DOWNSTREAM]);
        assertEquals(2, neData.Phases[FS_UPSTREAM]);
        assertEquals(1, neData.Phases[FS_DOWNSTREAM]);
        assertEquals(2, neData.ExonRank[FS_UPSTREAM]);
        assertEquals(4, neData.ExonRank[FS_DOWNSTREAM]);
        assertFalse(neData.phaseMatched());

        neData.setCodingBases(refGenome, 2);

        String upBases = chr1Bases.substring(25, 31);
        assertEquals(upBases, neData.CodingBases[FS_UPSTREAM]);

        String novelBases = chr1Bases.substring(40, 46); // + chr1Bases.substring(60, 71) + chr1Bases.substring(80, 86);
        assertEquals(novelBases, neData.NovelCodonBases);

        // String downBases = chr1Bases.substring(43, 51) + chr1Bases.substring(60, 61);
        assertEquals("", neData.CodingBases[FS_DOWNSTREAM]);

        // intronic to exonic - skips to next exon
        svData = new NeoEpitopeFusion(
                GENE_ID_1, GENE_ID_1, CHR_1, 35, POS_ORIENT, 1, GENE_ID_1, GENE_ID_1, CHR_1, 45, NEG_ORIENT,
                1, 1, "");

        neData = new SvNeoEpitope(svData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(CODING, neData.CodingType[FS_UPSTREAM]);
        assertEquals(CODING, neData.CodingType[FS_DOWNSTREAM]);
        assertEquals(INTRONIC, neData.RegionType[FS_UPSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWNSTREAM]);
        assertEquals(2, neData.Phases[FS_UPSTREAM]);
        assertEquals(1, neData.Phases[FS_DOWNSTREAM]);
        assertEquals(2, neData.ExonRank[FS_UPSTREAM]);
        assertEquals(4, neData.ExonRank[FS_DOWNSTREAM]);
        assertFalse(neData.phaseMatched());

        neData.setCodingBases(refGenome, 2);

        upBases = chr1Bases.substring(25, 31);
        assertEquals(upBases, neData.CodingBases[FS_UPSTREAM]);

        novelBases = chr1Bases.substring(40, 46); // + chr1Bases.substring(60, 71) + chr1Bases.substring(80, 86);
        assertEquals(novelBases, neData.NovelCodonBases);

        // String downBases = chr1Bases.substring(43, 51) + chr1Bases.substring(60, 61);
        assertEquals("", neData.CodingBases[FS_DOWNSTREAM]);

        // exonic to exonic, in phase
        svData = new NeoEpitopeFusion(
                GENE_ID_1, GENE_ID_1, CHR_1, 44, POS_ORIENT, 1, GENE_ID_1, GENE_ID_1, CHR_1, 63, NEG_ORIENT,
                1, 1, "");

        neData = new SvNeoEpitope(svData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(CODING, neData.CodingType[FS_UPSTREAM]);
        assertEquals(CODING, neData.CodingType[FS_DOWNSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_UPSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWNSTREAM]);
        assertEquals(1, neData.Phases[FS_UPSTREAM]);
        assertEquals(2, neData.Phases[FS_DOWNSTREAM]);
        assertEquals(3, neData.ExonRank[FS_UPSTREAM]);
        assertEquals(4, neData.ExonRank[FS_DOWNSTREAM]);
        assertTrue(neData.phaseMatched());

        neData.setCodingBases(refGenome, 2);

        upBases = chr1Bases.substring(28, 31) + chr1Bases.substring(40, 43);
        assertEquals(upBases, neData.CodingBases[FS_UPSTREAM]);

        novelBases = chr1Bases.substring(43, 45) + chr1Bases.substring(63, 64);
        assertEquals(novelBases, neData.NovelCodonBases);

        String downBases = chr1Bases.substring(64, 70);
        assertEquals(downBases, neData.CodingBases[FS_DOWNSTREAM]);

        // exonic to exonic with an insert sequence, out of phase
        String insSequence = "AA";

        svData = new NeoEpitopeFusion(
                GENE_ID_1, GENE_ID_1, CHR_1, 44, POS_ORIENT, 1, GENE_ID_1, GENE_ID_1, CHR_1, 63, NEG_ORIENT,
                1, 1, insSequence);

        neData = new SvNeoEpitope(svData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(CODING, neData.CodingType[FS_UPSTREAM]);
        assertEquals(CODING, neData.CodingType[FS_DOWNSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_UPSTREAM]);
        assertEquals(EXONIC, neData.RegionType[FS_DOWNSTREAM]);
        assertEquals(0, neData.Phases[FS_UPSTREAM]);
        assertEquals(2, neData.Phases[FS_DOWNSTREAM]);
        assertEquals(3, neData.ExonRank[FS_UPSTREAM]);
        assertEquals(4, neData.ExonRank[FS_DOWNSTREAM]);
        assertFalse(neData.phaseMatched());

        neData.setCodingBases(refGenome, 3);

        upBases = chr1Bases.substring(28, 31) + chr1Bases.substring(40, 45) + insSequence.substring(0, 1);
        assertEquals(upBases, neData.CodingBases[FS_UPSTREAM]);

        novelBases = insSequence.substring(1, 2) + chr1Bases.substring(63, 71) + chr1Bases.substring(80, 83);
        assertEquals(novelBases, neData.NovelCodonBases);

        assertEquals("", neData.CodingBases[FS_DOWNSTREAM]);

    }


    /*
    @Test
    public void testNeoEpitopes()
    {
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();

        tester.initialiseFusions(geneTransCache);

        PRE_GENE_PROMOTOR_DISTANCE = 10;

        // first a gene on the forward strand
        String geneName = "GENE1";
        String geneId1 = "ENSG0001";
        String chromosome1 = "1";
        byte strand = 1;

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId1, geneName, chromosome1, strand, 5, 100));
        addGeneData(geneTransCache, chromosome1, geneList);

        int transId = 1;

        int[] exonStarts = new int[]{5, 15, 25, 35, 45, 55, 65, 75};
        int[] exonPhases = new int[]{-1, 0, 2, 1, 0, 2, 1, -1};

        TranscriptData transData = createTransExons(geneId1, transId++, strand, exonStarts, exonPhases, 4, true);
        assertEquals(17, transData.CodingStart.longValue());

        List<TranscriptData> transDataList = Lists.newArrayList(transData);

        addTransExonData(geneTransCache, geneId1, transDataList);

        geneName = "GENE2";
        String geneId2 = "ENSG0002";
        String chromosome2 = "2";

        geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId2, geneName, chromosome2, strand, 5, 100));

        addGeneData(geneTransCache, chromosome2, geneList);

        transData = createTransExons(geneId2, transId++, strand, exonStarts, exonPhases, 4, true);
        transDataList = Lists.newArrayList(transData);

        addTransExonData(geneTransCache, geneId2, transDataList);

        MockRefGenome refGenome = new MockRefGenome();

        String refBases = "";
        String intron = "AAAAA";
        String nonCodingExon = "GGGGG";

        refBases += intron;
        refBases += nonCodingExon; // exon 1
        refBases += intron;
        refBases += "GGATG"; // exon 2, start of coding, ends on phase 2, so next is 0
        refBases += intron;
        refBases += "TCATC"; // exon 3, end phase 1
        refBases += intron;
        refBases += "ATCAT"; // exon 4, end phase 0
        refBases += intron;
        refBases += "CATCA"; // exon 5, end phase 2
        refBases += intron;
        refBases += "TCATC"; // exon 6, end phase 1
        refBases += intron;
        refBases += "ATCAT"; // exon 7, end phase 0
        refBases += intron;
        refBases += "CATAA"; // exon 8 including stop codon
        refBases += intron + intron;

        refGenome.RefGenomeMap.put(chromosome1, refBases);
        refGenome.RefGenomeMap.put(chromosome2, refBases);

        NeoEpitopeWriter neoEpFinder = new NeoEpitopeWriter(refGenome, geneTransCache, "");

        byte posOrient = 1;
        byte negOrient = -1;

        // create fusions between various phases sections of these transcripts
        // add upstream breakends
        List<GeneAnnotation> upGenes = Lists.newArrayList();
        int upPos = 22;
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, posOrient);

        // add downstream breakends
        List<GeneAnnotation> downGenes = Lists.newArrayList();
        int downPos = 22;
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, negOrient);

        FusionParameters params = new FusionParameters();

        // List<GeneFusion> fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params, false);
        List<GeneFusion> fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        NeoEpitopeFusion data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("M"));
        assertTrue(data.novelAcid().equals(""));
        assertTrue(data.downstreamAcids().equals("SSSSSSSSS"));

        // try again a phase 1 fusion
        upPos = 41;
        upGenes.clear();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, posOrient);

        // add downstream breakends
        downPos = 41;
        downGenes.clear();
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, negOrient);

        fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes,params);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("MSSS"));
        assertTrue(data.novelAcid().equals("S"));
        assertTrue(data.downstreamAcids().equals("SSSSS"));

        // phase 2 fusion
        upPos = 62;
        upGenes.clear();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, posOrient);

        // add downstream breakends
        downPos = 62;
        downGenes.clear();
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, negOrient);

        fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("MSSSSSS"));
        assertTrue(data.novelAcid().equals("S"));
        assertTrue(data.downstreamAcids().equals("SS"));

        // unphased fusion - only difference is that all downstream bases are collected (and may not form viable amino acids)
        upPos = 41;
        upGenes.clear();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, posOrient);

        // add downstream breakends
        downPos = 31;
        downGenes.clear();
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, negOrient);

        fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("MSSS"));
        assertTrue(data.novelAcid().equals("YHHHHHHH"));
        assertTrue(data.downstreamAcids().equals(""));
    }

    @Test
    public void testNeoEpitopesReverseStrand()
    {
        LinxTester tester = new LinxTester();

        EnsemblDataCache geneTransCache = createGeneDataCache();

        tester.initialiseFusions(geneTransCache);

        PRE_GENE_PROMOTOR_DISTANCE = 10;

        // first a gene on the forward strand
        String geneId1 = "ENSG0001";
        String chromosome1 = "1";
        byte strand = -1;

        int transId = 1;

        String intron = "AAAAA";
        String nonCodingExon = "GGGGG";
        String nonGenicDna = "TTTTT";
        String sCodon = convertAminoAcidToDnaCodon("S");

        // resultant bases and exon indices
        // AAAAACCCCCTTATGTTTTTATGATGATGATGATGATGATGATGATGATGATTTTTTGACATCCCCCAAAAA
        // 0    5    10   15   20                             51   56    62   67

        // first exon, starts coding
        String exon1 = nonCodingExon + swapRnaToDna(START_CODON) + sCodon;
        String transBases = exon1;

        // second exon - a string of 'S' amino acids following by 1 single base (ending on phase 1)
        String exon2 = "";
        int codonCount = 10;
        for(int i = 0; i < codonCount; ++i)
        {
            exon2 += sCodon;
        }

        // finish on phase 1
        exon2 += sCodon.substring(0, 1);

        transBases += intron + exon2;

        // final exon
        String exon3 = sCodon.substring(1) + swapRnaToDna(STOP_CODON_1) + nonCodingExon;
        transBases += intron + exon3;

        String revBases = reverseStrandBases(nonGenicDna + transBases + nonGenicDna);

        int transStart = nonGenicDna.length();
        int transEnd = transStart + transBases.length() - 1;

        // exons added from lower to higher positions
        List<ExonData> exons = Lists.newArrayList();
        exons.add(new ExonData(transId, transStart, transStart + exon3.length() - 1, 3, 1, -1));

        int nextExonStart = exons.get(exons.size() - 1).ExonEnd + intron.length() + 1;
        exons.add(new ExonData(transId, nextExonStart, nextExonStart + exon2.length() - 1, 2, 0, 1));

        nextExonStart = exons.get(exons.size() - 1).ExonEnd + intron.length() + 1;
        exons.add(new ExonData(transId, nextExonStart, nextExonStart + exon1.length() - 1, 1, -1, 0));

        TranscriptData transData = new TranscriptData(transId, generateTransName(transId), geneId1, true, strand, transStart, transEnd,
                transStart + nonCodingExon.length(), transEnd - nonCodingExon.length(), "");

        transData.exons().addAll(exons);

        addTransExonData(geneTransCache, geneId1, Lists.newArrayList(transData));

        // same format transcript on another chromosome, same strand
        String geneId2 = "ENSG0002";
        String chromosome2 = "2";

        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId1, "GENE1", chromosome1, strand, transStart, transEnd));
        addGeneData(geneTransCache, chromosome1, geneList);

        geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId2, "GENE2", chromosome2, strand, transStart, transEnd));
        addGeneData(geneTransCache, chromosome2, geneList);

        transData = new TranscriptData(++transId, generateTransName(transId), geneId2, true, strand, transStart, transEnd,
                transStart + nonCodingExon.length(), transEnd - nonCodingExon.length(), "");

        transData.exons().addAll(exons);

        addTransExonData(geneTransCache, geneId2, Lists.newArrayList(transData));

        MockRefGenome refGenome = new MockRefGenome();

        refGenome.RefGenomeMap.put(chromosome1, revBases);
        refGenome.RefGenomeMap.put(chromosome2, revBases);

        NeoEpitopeWriter neoEpFinder = new NeoEpitopeWriter(refGenome, geneTransCache, "");

        byte posOrient = 1;
        byte negOrient = -1;

        // create fusions between various phases sections of these transcripts
        // add upstream breakends
        List<GeneAnnotation> upGenes = Lists.newArrayList();
        int upPos = 51;
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, posOrient);

        // add downstream breakends
        List<GeneAnnotation> downGenes = Lists.newArrayList();
        int downPos = 51;
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, negOrient);

        FusionParameters params = new FusionParameters();

        List<GeneFusion> fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        NeoEpitopeFusion data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("MS"));
        assertTrue(data.novelAcid().equals(""));
        assertTrue(data.downstreamAcids().equals("SSSSSSSSSS"));

        // with a fusion between the 2nd and 3rd exons, splitting a codon
        upPos = 17;
        upGenes.clear();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, posOrient);

        // add downstream breakends
        downGenes.clear();
        downPos = 17;
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, negOrient);

        fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("SSSSSSSSSS"));
        assertTrue(data.novelAcid().equals("S"));
        assertTrue(data.downstreamAcids().equals("_"));

        // a fusion to the 5'UTR of the downstream gene
        upPos = 51;
        upGenes.clear();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, negOrient);

        // add downstream breakends
        downGenes.clear();
        downPos = 68;
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, posOrient);

        fusions = tester.FusionAnalyser.getFusionFinder().findFusions(upGenes, downGenes, params);

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("MS"));
        assertTrue(data.novelAcid().equals(""));
        assertTrue(data.downstreamAcids().equals("SSSSSSSSSS"));

        // a fusion to the 3'UTR of the downstream gene
        upPos = 51;
        upGenes.clear();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, negOrient);

        // add downstream breakends
        downGenes.clear();
        downPos = 14;
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome2, downPos, posOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome2, downPos, posOrient);

        fusions.clear();

        neoEpFinder.checkFusions(fusions, upGenes, downGenes);
        assertEquals(1, fusions.size());

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("MS"));
        assertTrue(data.novelAcid().equals("H"));
        assertTrue(data.downstreamAcids().equals(""));

        // again with a longer downstream 3'UTR sequence
        String geneId3 = "ENSG0003";
        String chromosome3 = "3";
        strand = 1;

        // resultant bases and exon indices
        // AAAAACCCCCTTATGTTTTTATGATGATGATGATGATGATGATGATGATGATTTTTTGACATCCCCCAAAAA
        // 0    5    10   15   20                             51   56    62   67

        transBases = nonCodingExon; // exon 1
        transBases += intron + nonCodingExon + swapRnaToDna(START_CODON) + sCodon; // exon 2
        transBases += intron + sCodon + swapRnaToDna(STOP_CODON_1) + nonCodingExon; // exon 3
        transBases += intron + nonCodingExon; // exon 4
        transBases += intron + nonCodingExon; // exon 5
        transBases += intron + nonCodingExon; // exon 6

        String refBases = nonGenicDna + transBases + nonGenicDna;

        refGenome.RefGenomeMap.put(chromosome3, refBases);

        transStart = nonGenicDna.length();
        transEnd = transStart + transBases.length() - 1;


        exons.clear();
        exons.add(new ExonData(transId, transStart, transStart + nonCodingExon.length() - 1, 1, -1, -1));

        nextExonStart = exons.get(exons.size() - 1).ExonEnd + intron.length() + 1;
        exons.add(new ExonData(transId, nextExonStart, nextExonStart + nonCodingExon.length() + 6 - 1, 2, -1, 0));

        nextExonStart = exons.get(exons.size() - 1).ExonEnd + intron.length() + 1;
        exons.add(new ExonData(transId, nextExonStart, nextExonStart + 6 + nonCodingExon.length() - 1, 3, 0, -1));

        nextExonStart = exons.get(exons.size() - 1).ExonEnd + intron.length() + 1;
        exons.add(new ExonData(transId, nextExonStart, nextExonStart + nonCodingExon.length() - 1, 4, -1, -1));

        nextExonStart = exons.get(exons.size() - 1).ExonEnd + intron.length() + 1;
        exons.add(new ExonData(transId, nextExonStart, nextExonStart + nonCodingExon.length() - 1, 5, -1, -1));

        nextExonStart = exons.get(exons.size() - 1).ExonEnd + intron.length() + 1;
        exons.add(new ExonData(transId, nextExonStart, nextExonStart + nonCodingExon.length() - 1, 6, -1, -1));

        int codingStart = exons.get(1).ExonStart + nonCodingExon.length() - 1;
        int codingEnd = exons.get(2).ExonStart + 6 - 1;

        transData = new TranscriptData(++transId, generateTransName(transId), geneId3, true, strand, transStart, transEnd,
                codingStart, codingEnd, "");

        transData.exons().addAll(exons);

        geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(geneId3, "GENE3", chromosome3, strand, transStart, transEnd));
        addGeneData(geneTransCache, chromosome3, geneList);

        addTransExonData(geneTransCache, geneId3, Lists.newArrayList(transData));

        // a fusion to the 3'UTR of the downstream gene
        upPos = 51;
        upGenes.clear();
        upGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, true, chromosome1, upPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        upGenes.get(0).setPositionalData(chromosome1, upPos, negOrient);

        // add downstream breakends
        downGenes.clear();
        downPos = codingEnd + nonCodingExon.length() + 2;
        downGenes.addAll(geneTransCache.findGeneAnnotationsBySv(0, false, chromosome3, downPos, negOrient, PRE_GENE_PROMOTOR_DISTANCE));
        downGenes.get(0).setPositionalData(chromosome3, downPos, negOrient);

        fusions.clear();

        neoEpFinder.checkFusions(fusions, upGenes, downGenes);
        assertEquals(1, fusions.size());

        neoEpFinder.reportNeoEpitopes(tester.SampleId, fusions);

        assertEquals(1, neoEpFinder.getResults().size());
        data = neoEpFinder.getResults().get(0);
        assertTrue(data.upstreamAcids().equals("MS"));
        assertTrue(data.novelAcid().equals("GGGGG"));
        assertTrue(data.downstreamAcids().equals(""));
        assertEquals(10, data.downstreamNmdBases());
    }
    */

}
