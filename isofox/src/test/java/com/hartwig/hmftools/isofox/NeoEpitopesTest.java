package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.fusionInfo;
import static com.hartwig.hmftools.common.neo.NeoEpitopeType.INFRAME_FUSION;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_1;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_2;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_2;
import static com.hartwig.hmftools.isofox.TestUtils.createCigar;
import static com.hartwig.hmftools.isofox.TestUtils.createReadRecord;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentMatcher.calcBaseOverlap;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentMatcher.calcCoordinatesOverlap;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentMatcher.expandRange;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentMatcher.findFusionSupport;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentSupport.EXACT_MATCH;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;

import java.util.List;

import com.hartwig.hmftools.common.neo.NeoEpitopeFile;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.neo.NeoEpitopeData;
import com.hartwig.hmftools.isofox.neo.NeoFragmentSupport;

import org.apache.commons.compress.utils.Lists;
import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class NeoEpitopesTest
{
    @Test
    public void testNovelRanges()
    {
        List<int[]> coords = Lists.newArrayList();

        coords.add(new int[] {10, 20});
        coords.add(new int[] {30, 40});
        coords.add(new int[] {50, 60});
        coords.add(new int[] {70, 80});

        List<int[]> novelRanges = Lists.newArrayList();
        novelRanges.add(new int[] { 35, 35 });
        expandRange(novelRanges, 35, coords, 5, false);
        expandRange(novelRanges, 35, coords, 5, true);

        assertEquals(1, novelRanges.size());
        assertEquals(30, novelRanges.get(0)[SE_START]);
        assertEquals(40, novelRanges.get(0)[SE_END]);

        expandRange(novelRanges, 30, coords, 10, false);
        expandRange(novelRanges, 40, coords, 17, true);

        assertEquals(4, novelRanges.size());
        assertEquals(11, novelRanges.get(0)[SE_START]);
        assertEquals(75, novelRanges.get(3)[SE_END]);
    }

    @Test
    public void testCoordsOverlap()
    {

        int[] range1 = new int[] {10, 20};
        int[] range2 = new int[] {30, 40};
        assertEquals(0, calcBaseOverlap(range1, range2));

        range2 = new int[] {20, 30};
        assertEquals(1, calcBaseOverlap(range1, range2));

        range2 = new int[] {0, 11};
        assertEquals(2, calcBaseOverlap(range1, range2));

        range2 = new int[] {5, 25};
        assertEquals(11, calcBaseOverlap(range1, range2));

        List<int[]> coords1 = Lists.newArrayList();
        coords1.add(new int[] {10, 20});
        coords1.add(new int[] {30, 40});
        coords1.add(new int[] {50, 60});
        coords1.add(new int[] {70, 80});

        assertEquals(44, calcCoordinatesOverlap(coords1, coords1));

        List<int[]> coords2 = Lists.newArrayList();
        coords2.add(new int[] {0, 10});
        coords2.add(new int[] {21, 30});
        coords2.add(new int[] {41, 50});

        assertEquals(3, calcCoordinatesOverlap(coords1, coords2));
    }

    @Test
    public void testFusions()
    {
        String codingBases = generateRandomBases(310);

        // fusion DEL from 50 to 250 bases
        int cbUpStart = 20;
        int cbUpEnd = 50;
        int cbDownStart = 250;
        int cbDownEnd = 280;

        NeoEpitopeData neData = createFusionNeoEpitope(
                INFRAME_FUSION, CHR_1, cbUpEnd, POS_ORIENT, CHR_1, cbDownStart, NEG_ORIENT,
                GENE_ID_1, GENE_ID_2, "", "",
                cbUpStart, cbUpEnd, codingBases.substring(cbUpStart, cbUpEnd), createCigar(0, cbUpEnd - cbUpStart + 1, 0).toString(),
                cbDownStart, cbDownEnd, codingBases.substring(cbDownStart, cbDownEnd), createCigar(0,
                        cbDownEnd - cbDownStart + 1, 0).toString());

        // first read has no overlap
        int rbStart = 80;
        int rbEnd = 100;
        ReadRecord read = createReadRecord(1, CHR_1, rbStart, rbEnd, codingBases.substring(rbStart, rbEnd + 1),
                createCigar(0, rbEnd - rbStart + 1, 0));

        NeoFragmentSupport support = findFusionSupport(neData, FS_UP, read);

        assertFalse(support.hasAnySupport());

        // upstream exact-match overlap
        rbStart = 10;
        rbEnd = 40;
        read = createReadRecord(1, CHR_1, rbStart, rbEnd, codingBases.substring(rbStart, rbEnd + 1),
                createCigar(0, rbEnd - rbStart + 1, 0));

        support = findFusionSupport(neData, FS_UP, read);

        assertEquals(1, support.UpFragments[EXACT_MATCH]);

        // downstream exact-match overlap
        rbStart = 260;
        rbEnd = 290;
        read = createReadRecord(1, CHR_1, rbStart, rbEnd, codingBases.substring(rbStart, rbEnd + 1),
                createCigar(0, rbEnd - rbStart + 1, 0));

        support = findFusionSupport(neData, FS_DOWN, read);

        assertEquals(1, support.DownFragments[EXACT_MATCH]);

        // novel support from upstream side
        rbStart = 30;
        rbEnd = 50;
        read = createReadRecord(1, CHR_1, rbStart, rbEnd,
                codingBases.substring(rbStart, rbEnd + 1) + codingBases.substring(250, 260),
                createCigar(0, rbEnd - rbStart + 1, 10));

        support = findFusionSupport(neData, FS_UP, read);

        assertEquals(1, support.NovelFragments[EXACT_MATCH]);

        // novel support from downstream side
        rbStart = 250;
        rbEnd = 270;
        read = createReadRecord(1, CHR_1, rbStart, rbEnd,
                codingBases.substring(40, 50) + codingBases.substring(rbStart, rbEnd + 1),
                createCigar(10, rbEnd - rbStart + 1, 0));

        support = findFusionSupport(neData, FS_DOWN, read);

        assertEquals(1, support.NovelFragments[EXACT_MATCH]);

        // neg-strand to neg strand fusion on different chromosome, with 19 base N-splits in each
        // upstream intron - 261-279, downstream intron - 21-39
        cbUpStart = 250;
        cbUpEnd = 290;
        cbDownStart = 10;
        cbDownEnd = 50;

        neData = createFusionNeoEpitope(
                INFRAME_FUSION, CHR_1, cbUpStart, NEG_ORIENT, CHR_2, cbDownStart, NEG_ORIENT,
                GENE_ID_1, GENE_ID_2, "", "",
                cbUpStart, cbUpEnd, codingBases.substring(250, 261) + codingBases.substring(280, 291),
                createCigar(0, 11, 19, 11, 0).toString(),
                cbDownStart, cbDownEnd, codingBases.substring(10, 21) + codingBases.substring(40, 51),
                createCigar(0, 11, 19, 11, 0).toString());

        // upstream exact-match overlap
        rbStart = 255;
        rbEnd = 295;
        read = createReadRecord(1, CHR_1, rbStart, rbEnd,
                codingBases.substring(255, 261) + codingBases.substring(280, 297),
                createCigar(0, 6, 19, 17, 0));

        support = findFusionSupport(neData, FS_UP, read);

        assertEquals(1, support.UpFragments[EXACT_MATCH]);

        // invalid upstream exact-match overlap
        rbStart = 245;
        rbEnd = 295;
        read = createReadRecord(1, CHR_1, rbStart, rbEnd,
                codingBases.substring(245, 261) + codingBases.substring(280, 297),
                createCigar(0, 11, 19, 17, 0));

        support = findFusionSupport(neData, FS_UP, read);

        assertFalse(support.hasAnySupport());

        // downstream exact-match overlap
        rbStart = 15;
        rbEnd = 60;
        read = createReadRecord(1, CHR_1, rbStart, rbEnd,
                codingBases.substring(15, 21) + codingBases.substring(40, 61),
                createCigar(0, 6, 19, 21, 0));

        support = findFusionSupport(neData, FS_DOWN, read);

        assertEquals(1, support.DownFragments[EXACT_MATCH]);

        // invalid downstream exact-match overlap
        rbStart = 5; // starts before the junction
        rbEnd = 60;
        read = createReadRecord(1, CHR_1, rbStart, rbEnd,
                codingBases.substring(5, 21) + codingBases.substring(40, 61),
                createCigar(0, 6, 19, 21, 0));

        support = findFusionSupport(neData, FS_DOWN, read);

        assertFalse(support.hasAnySupport());

        // read crossing the junction matching upstream read bases, downstream soft-clipped
        rbStart = 250;
        rbEnd = 290;
        read = createReadRecord(1, CHR_1, rbStart, rbEnd,
                reverseComplementBases(codingBases.substring(10, 19)) + codingBases.substring(250, 261) + codingBases.substring(280, 291),
                createCigar(10, 11, 19, 11, 0));

        support = findFusionSupport(neData, FS_UP, read);

        // read crossing the junction matching downstream read bases, upstream soft-clipped
        assertEquals(1, support.NovelFragments[EXACT_MATCH]);

        rbStart = 10;
        rbEnd = 50;
        read = createReadRecord(1, CHR_1, rbStart, rbEnd,
                reverseComplementBases(codingBases.substring(250, 260)) + codingBases.substring(10, 21) + codingBases.substring(40, 51),
                createCigar(10, 11, 19, 11, 0));

        support = findFusionSupport(neData, FS_DOWN, read);

        assertEquals(1, support.NovelFragments[EXACT_MATCH]);

        // N-split fusion on the neg-strand
        cbUpStart = 250;
        cbUpEnd = 290;
        cbDownStart = 10;
        cbDownEnd = 50;

        neData = createFusionNeoEpitope(
                INFRAME_FUSION, CHR_1, cbUpStart, NEG_ORIENT, CHR_1, cbDownEnd, POS_ORIENT,
                GENE_ID_1, GENE_ID_2, "", "",
                cbUpStart, cbUpEnd, codingBases.substring(250, 261) + codingBases.substring(280, 291),
                createCigar(0, 11, 19, 11, 0).toString(),
                cbDownStart, cbDownEnd, codingBases.substring(10, 21) + codingBases.substring(40, 51),
                createCigar(0, 11, 19, 11, 0).toString());

        // read crossing the junction matching upstream read bases, downstream soft-clipped
        rbStart = 10;
        rbEnd = 290;

        Cigar readCigar = new Cigar();
        readCigar.add(new CigarElement(20 - 10 + 1, CigarOperator.M));
        readCigar.add(new CigarElement(39 - 21 + 1, CigarOperator.N));
        readCigar.add(new CigarElement(50 - 40 + 1, CigarOperator.M));
        readCigar.add(new CigarElement(249 - 51 + 1, CigarOperator.N)); // fusion gap
        readCigar.add(new CigarElement(260 - 250 + 1, CigarOperator.M));
        readCigar.add(new CigarElement(279 - 261, CigarOperator.N));
        readCigar.add(new CigarElement(290 - 280 + 1, CigarOperator.M));

        read = createReadRecord(1, CHR_1, rbStart, rbEnd,
                codingBases.substring(10, 21) + codingBases.substring(40, 51)
                        + codingBases.substring(250, 261) + codingBases.substring(280, 291),
                readCigar);

        support = findFusionSupport(neData, FS_UP, read);
        assertEquals(1, support.NovelFragments[EXACT_MATCH]);
    }

    private NeoEpitopeData createFusionNeoEpitope(
            final NeoEpitopeType varType, final String chrUp, int posUp, byte orientUp, final String chrDown, int posDown, byte orientDown,
            final String geneIdUp, final String geneIdDown, final String transcriptsUp, final String transcriptsDown,
            int codingBaseUpPosStart, int codingBaseUpPosEnd, final String codingBasesUp, final String codingBaseCigarUp,
            int codingBaseDownPosStart, int codingBaseDownPosEnd, final String codingBasesDown, final String codingBaseCigarDown)
    {
        final String varInfo = fusionInfo(new String[] {chrUp, chrDown}, new int[] { posUp, posDown}, new byte[] { orientUp, orientDown });

        final NeoEpitopeFile neFile = new NeoEpitopeFile(
                1, varType, varInfo, 1, 1, 0, geneIdUp, geneIdDown, geneIdUp, geneIdDown,
                chrUp, chrDown, orientUp, orientDown,
                "", "", "", 0, 0, 0, 0,
                0, 0, 0,
                transcriptsUp, transcriptsDown, "",
                codingBaseUpPosStart, codingBaseUpPosEnd, codingBasesUp, codingBaseCigarUp,
                codingBaseDownPosStart, codingBaseDownPosEnd, codingBasesDown, codingBaseCigarDown);

        return new NeoEpitopeData(neFile);
    }
}
