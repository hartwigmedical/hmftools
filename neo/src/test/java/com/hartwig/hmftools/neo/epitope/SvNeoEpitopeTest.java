package com.hartwig.hmftools.neo.epitope;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.neo.epitope.SvNeoEpitope.svIsNonDisruptiveInCodingTranscript;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertFalse;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Assert;
import org.junit.Test;

public class SvNeoEpitopeTest
{
    @Test
    public void testSvFusions()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        final String chr1Bases = generateRandomBases(120);

        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        // tests: intronic vs exonic, with and without insertions

        int[] exonStarts = { 0, 20, 40, 60, 80, 100 };
        Integer codingStart = Integer.valueOf(25);
        Integer codingEnd = Integer.valueOf(85);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        final String[] validTrans = new String[] {transDataPosStrand.TransName, transDataPosStrand.TransName};

        // intronic
        NeoEpitopeFusion svData = new NeoEpitopeFusion(
                GENE_ID_1, GENE_ID_1, CHR_1, 35, ORIENT_FWD, 1, GENE_ID_1, GENE_ID_1, CHR_1, 55, ORIENT_REV,
                1, 1, 1, "", 0, validTrans);

        NeoEpitope neData = new SvNeoEpitope(svData);

        neData.setTranscriptData(transDataPosStrand, transDataPosStrand);

        assertEquals(35, neData.position(FS_UP));
        assertEquals(55, neData.position(FS_DOWN));
        assertEquals(ORIENT_FWD, neData.orientation(FS_UP));
        assertEquals(ORIENT_REV, neData.orientation(FS_DOWN));

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
        String upWildAAs = EpitopeUtils.getAminoAcids(upWildtypeBases, true);
        assertEquals(upWildAAs, neData.UpstreamWildTypeAcids);

        // intronic to exonic - skips to next exon
        svData = new NeoEpitopeFusion(
                GENE_ID_1, GENE_ID_1, CHR_1, 35, ORIENT_FWD, 1, GENE_ID_1, GENE_ID_1, CHR_1, 45, ORIENT_REV,
                1, 1, 1, "", 0, validTrans);

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
                GENE_ID_1, GENE_ID_1, CHR_1, 44, ORIENT_FWD, 1, GENE_ID_1, GENE_ID_1, CHR_1, 63, ORIENT_REV,
                1, 1, 1, "", 0, validTrans);

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
        assertTrue(neData.phaseMatched());

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
                GENE_ID_1, GENE_ID_1, CHR_1, 44, ORIENT_FWD, 1, GENE_ID_1, GENE_ID_1, CHR_1, 63, ORIENT_REV,
                1, 1,1, insSequence, 0, validTrans);

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
        upWildAAs = EpitopeUtils.getAminoAcids(upWildtypeBases, true);
        assertEquals(upWildAAs, neData.UpstreamWildTypeAcids);

        // intronic to upstream, skips to first splice acceptor

    }

    @Test
    public void testNonDisruptiveSv()
    {
        int[] exonStarts = { 100, 300, 500, 700, 900 };
        Integer codingStart = Integer.valueOf(150);
        Integer codingEnd = Integer.valueOf(750);

        TranscriptData transData1 = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 100, codingStart, codingEnd, false, "");

        final String[] validTrans = new String[] {transData1.TransName, transData1.TransName};

        // intronic DEL
        int[] svPositions = { 420, 480 };

        assertTrue(svIsNonDisruptiveInCodingTranscript(svPositions, transData1));

        svPositions = new int[] { 320, 480 }; // exonic
        assertFalse(svIsNonDisruptiveInCodingTranscript(svPositions, transData1));

        svPositions = new int[] { 420, 680 }; // deletes an exon
        assertFalse(svIsNonDisruptiveInCodingTranscript(svPositions, transData1));

        svPositions = new int[] { 820, 880 }; // UTR region
        assertFalse(svIsNonDisruptiveInCodingTranscript(svPositions, transData1));
    }

}
