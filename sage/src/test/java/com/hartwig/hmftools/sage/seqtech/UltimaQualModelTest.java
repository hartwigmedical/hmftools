package com.hartwig.hmftools.sage.seqtech;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.common.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.sage.common.TestUtils.setIlluminaSequencing;
import static com.hartwig.hmftools.sage.common.TestUtils.setUltimaSequencing;
import static com.hartwig.hmftools.sage.seqtech.UltimaModelType.BASE_SHIFT;
import static com.hartwig.hmftools.sage.seqtech.UltimaModelType.HOMOPOLYMER_ADJUSTMENT;
import static com.hartwig.hmftools.sage.seqtech.UltimaModelType.HOMOPOLYMER_DELETION;
import static com.hartwig.hmftools.sage.seqtech.UltimaModelType.HOMOPOLYMER_TRANSITION;
import static com.hartwig.hmftools.sage.seqtech.UltimaModelType.OTHER;
import static com.hartwig.hmftools.sage.seqtech.UltimaModelType.SNV;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.BQR_CACHE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;

import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SeqTechTestUtils;
import com.hartwig.hmftools.common.variant.SimpleVariant;

import org.junit.After;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class UltimaQualModelTest
{
    private final MockRefGenome mRefGenome;
    private final UltimaQualModelBuilder mModelBuilder;

    private static final String BUFFER_REF_BASES = "AACCGGTTAACCGGTTAA";

    public UltimaQualModelTest()
    {
        mRefGenome = new MockRefGenome(true);
        mModelBuilder = new UltimaQualModelBuilder(mRefGenome);
        setUltimaSequencing();
        populateUltimaBqrCache();
    }

    protected static void populateUltimaBqrCache()
    {
        // ultima_bqr_recalibration.tsv
        // InputStream inputStream = UltimaQualModelTest.class.getResourceAsStream("/seqtech/ultima_bqr_recalibration.tsv");
        // List<String> lines = new BufferedReader(new InputStreamReader(inputStream)).lines().collect(toList());
        // BQR_CACHE.loadRecalibrationData(lines);

        UltimaUtils.loadBqrCache(null);
    }

    @After
    public void resetSequencingType() { setIlluminaSequencing(); }

    private void setRefBases(final String refBases)
    {
        mRefGenome.RefGenomeMap.put(CHR_1, refBases);
    }

    @Test
    public void testHomopolymerAdjustment()
    {
        //                                    0123456789
        String refBases = BUFFER_REF_BASES + "ATTTTCGTCGT" + BUFFER_REF_BASES;
        setRefBases(refBases);
        // delete of 1 base
        SimpleVariant variant = new SimpleVariant(CHR_1, 19, "AT", "A");
        UltimaQualModel model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertNotNull(model);
        assertEquals(HOMOPOLYMER_ADJUSTMENT, model.type());

        String readBases = BUFFER_REF_BASES + "ATTTCGTCGT";
        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());
        byte[] t0Values = buildDefaultBaseQuals(readBases.length());
        byte[] tpValues = new byte[readBases.length()];
        tpValues[19] = 1;
        tpValues[21] = 1;
        SAMRecord read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        byte calcQual = model.calculateQual(read, 18);
        assertEquals(38, calcQual);

        tpValues = new byte[readBases.length()];
        tpValues[20] = 1;
        baseQualities[20] = 25;
        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 18);
        assertEquals(25, calcQual);

        // delete of 3 bases
        variant = new SimpleVariant(CHR_1, 19, "ATTT", "A");

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertEquals(HOMOPOLYMER_ADJUSTMENT, model.type());

        readBases = BUFFER_REF_BASES + "ATCGTCGT";
        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];
        tpValues[19] = 3;
        baseQualities[19] = 25;
        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 18);
        assertEquals(25, calcQual);

        // test where the required adjustment value isn't in the TP values
        tpValues[19] = -1;

        calcQual = model.calculateQual(read, 18);

        // ULTIMA TODO
        // assertEquals(Math.min(49, ULTIMA_MAX_QUAL_TP), calcQual);

        // insert of 2 bases
        variant = new SimpleVariant(CHR_1, 19, "A", "ATT");

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertEquals(HOMOPOLYMER_ADJUSTMENT, model.type());

        readBases = BUFFER_REF_BASES + "ATTTTTTCGTCGT";
        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];
        tpValues[19] = 1;
        tpValues[20] = -1;
        tpValues[21] = -2; // will use these 2
        tpValues[22] = -2;
        baseQualities[21] = 30;
        baseQualities[22] = 30;
        tpValues[23] = 11;
        tpValues[24] = 1;
        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 18);
        assertEquals(27, calcQual);

        // C>CA in TT C GTC, insert of base which doesn't match the existing context
        variant = new SimpleVariant(CHR_1, 24, "C", "CA");

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertEquals(HOMOPOLYMER_ADJUSTMENT, model.type());

        readBases = BUFFER_REF_BASES + "ATTTTCAGTCGT" + BUFFER_REF_BASES;

        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];

        tpValues[25] = -1;
        baseQualities[25] = 32;

        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 24);
        assertEquals(32, calcQual);
    }

    @Test
    public void testHomopolymerDeletion()
    {
        //                                    0123456789
        String refBases = BUFFER_REF_BASES + "ATTTTCGACGT" + BUFFER_REF_BASES;
        setRefBases(refBases);

        // the whole HP must be deleted
        SimpleVariant variant = new SimpleVariant(CHR_1, 19, "ATTTT", "A");

        UltimaQualModel model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertNotNull(model);
        assertEquals(HOMOPOLYMER_DELETION, model.type());

        String readBases = BUFFER_REF_BASES + "ACGTCGT";
        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());
        byte[] t0Values = buildDefaultBaseQuals(readBases.length());
        byte[] tpValues = new byte[readBases.length()];
        t0Values[18] = 25;
        t0Values[19] = 30;
        SAMRecord read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        byte calcQual = model.calculateQual(read, 18);
        assertEquals(30, calcQual);

        read.setReadNegativeStrandFlag(true);

        calcQual = model.calculateQual(read, 18);
        assertEquals(30, calcQual);

        // test out of cycle deletions
        variant = new SimpleVariant(CHR_1, 25, "GA", "G");

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertNotNull(model);
        assertEquals(HOMOPOLYMER_DELETION, model.type());

        readBases = BUFFER_REF_BASES + "ATTTTCGCGT";
        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];
        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 18);
        assertEquals(BQR_CACHE.outOfCycleT0Qual(), calcQual);

        read.setReadNegativeStrandFlag(true);
        calcQual = model.calculateQual(read, 18);
        assertEquals(BQR_CACHE.outOfCycleT0Qual(), calcQual);
    }

    @Test
    public void testHomopolymerTransition()
    {
        //                                    01234567890
        String refBases = BUFFER_REF_BASES + "ATTTTAAAAAG" + BUFFER_REF_BASES;
        setRefBases(refBases);

        // a deletion which crosses 2 HPs
        SimpleVariant variant = new SimpleVariant(CHR_1, 21, "TTTAA", "T");

        UltimaQualModel model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertNotNull(model);
        assertEquals(HOMOPOLYMER_TRANSITION, model.type());

        //                                     0123456
        String readBases = BUFFER_REF_BASES + "ATTAAAG";
        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());
        byte[] t0Values = buildDefaultBaseQuals(readBases.length());
        byte[] tpValues = new byte[readBases.length()];

        tpValues[20] = 2;
        tpValues[21] = 2;
        baseQualities[20] = 11;
        baseQualities[21] = 11;

        tpValues[22] = 2;
        tpValues[24] = 2;
        baseQualities[22] = 25;
        baseQualities[24] = 25;

        SAMRecord read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        byte calcQual = model.calculateQual(read, 21);
        assertEquals(30, calcQual);

        // as before but with longer, lopsided transitional delete
        variant = new SimpleVariant(CHR_1, 20, "TTTTA", "T");

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertEquals(HOMOPOLYMER_TRANSITION, model.type());

        //                              0123456
        readBases = BUFFER_REF_BASES + "ATAAAAG";
        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];

        tpValues[20] = 3;
        baseQualities[20] = 25;

        tpValues[21] = -1;
        tpValues[22] = 1;
        tpValues[23] = 1;
        tpValues[24] = -1;
        baseQualities[22] = 12;
        baseQualities[23] = 12;

        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 20);
        assertEquals(34, calcQual);
    }

    @Test
    public void testSNVs()
    {
        //                                    01     234     56
        String refBases = BUFFER_REF_BASES + "AG" + "ACA" + "AG" + BUFFER_REF_BASES;
        setRefBases(refBases);

        // C>T in ACA > ATA, matches neither side's ref so reverts to max qual
        SimpleVariant variant = new SimpleVariant(CHR_1, 22, "C", "T");

        UltimaQualModel model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertNotNull(model);
        assertEquals(SNV, model.type());

        String readBases = BUFFER_REF_BASES + "AGATAAG";
        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());
        byte[] t0Values = buildDefaultBaseQuals(readBases.length());
        byte[] tpValues = new byte[readBases.length()];

        SAMRecord read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        byte calcQual = model.calculateQual(read, 22);
        assertEquals(BQR_CACHE.outOfCycleT0Qual(), calcQual);

        //                             01     234     56
        refBases = BUFFER_REF_BASES + "AG" + "CCG" + "AG" + BUFFER_REF_BASES;
        setRefBases(refBases);

        // C>T in CCG > CTG, left contraction, right insertion/expansion
        variant = new SimpleVariant(CHR_1, 22, "C", "T");

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));

        readBases = BUFFER_REF_BASES + "AGCTGAG";

        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];

        // the deleted base qual
        tpValues[21] = 1;
        baseQualities[21] = 16;

        tpValues[22] = -1;
        baseQualities[22] = 21;

        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 22);
        assertEquals(21, calcQual);

        // C>T GCT > GTT, left full delete, right insertion/expansion
        //                             01     234     56
        refBases = BUFFER_REF_BASES + "AG" + "GCT" + "AG" + BUFFER_REF_BASES;
        setRefBases(refBases);

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));

        //                              0123456
        readBases = BUFFER_REF_BASES + "AGGTTAG" + BUFFER_REF_BASES;

        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];

        t0Values[21] = 10;
        t0Values[22] = 15; // an SNV with any HP deletion, including of a single base, just uses t0

        tpValues[22] = -1;
        baseQualities[22] = 21;

        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 22);
        assertEquals(21, calcQual);

        // C>T in TCA > TTA, left 1-base ins/expansion, right full delete
        //                             01     234     56
        refBases = BUFFER_REF_BASES + "AG" + "TCA" + "AG" + BUFFER_REF_BASES;
        setRefBases(refBases);

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));

        //                              0123456
        readBases = BUFFER_REF_BASES + "AGTTAAG" + BUFFER_REF_BASES;

        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];

        t0Values[22] = 25;
        t0Values[23] = 10;

        tpValues[21] = tpValues[22] = -1;
        baseQualities[21] = baseQualities[22] = 20;

        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 22);
        assertEquals(25, calcQual);


        // C>T in GCC > GTC, left contraction, right insertion/expansion
        //                             01     234     56
        refBases = BUFFER_REF_BASES + "AG" + "GCC" + "AG" + BUFFER_REF_BASES;
        setRefBases(refBases);

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));

        readBases = BUFFER_REF_BASES + "AGGTCAG";

        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];

        // the delete base qual
        tpValues[23] = 1;
        baseQualities[23] = 16;

        tpValues[22] = -1;
        baseQualities[22] = 21;

        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 22);
        assertEquals(21, calcQual);


        // contraction on the C side, expansion on the T side
        // C>T in TTT TCC CCC > TTT TTC CCC

        variant = new SimpleVariant(CHR_1, 23, "C", "T");

        //                             012     345     678
        refBases = BUFFER_REF_BASES + "TTT" + "TCC" + "CCC" + BUFFER_REF_BASES;
        setRefBases(refBases);

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));

        //                              012345678
        readBases = BUFFER_REF_BASES + "TTTTTCCCC" + BUFFER_REF_BASES;

        baseQualities = buildDefaultBaseQuals(readBases.length());
        t0Values = buildDefaultBaseQuals(readBases.length());
        tpValues = new byte[readBases.length()];

        // the HP insert base quals
        tpValues[19] = tpValues[23] = -1;
        baseQualities[19] = baseQualities[23] = 16;

        // the HP delete quals
        tpValues[24] = tpValues[27] = 1;
        baseQualities[24] = baseQualities[27] = 21;

        read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        calcQual = model.calculateQual(read, 23);
        assertEquals(18, calcQual);

        // HP on right partially deleted -> right matches ref, insert on left
        // AG GCC AGA > AG GTC AGA

        // HP on at var position fully deleted -> right matches alt, like a HP expansion on right
        // AG GCT AGA > AG GTT AGA

        // HP at var position fully deleted -> left matches alt, like a HP expansion on left
        // AG TCG AGA > AG TTG AGA


        // variant has one straddling base matching ref or alt
        // AG CCG AGA > AG CTG AGA
        // delete of C on left,  insert of T


        // variant has one straddling base matching ref or alt
        // TC TCG GC > TC TTG GC
        // full delete of C, insert of additional T on left
    }

    @Test
    public void testBaseShifts()
    {
        // position starting at 19: 9012 345 6789
        // index                    0134 567 8901
        // ref bases                CAAA ACC TTTC
        // read bases               CAAA CCG TTTC
        String leftRefBases = "CAAA";
        String varRefBases = "ACC";
        String varAltBases = "CCG";
        String rightRefBases = "TTTC";
        String refBases = BUFFER_REF_BASES + leftRefBases + varRefBases + rightRefBases + BUFFER_REF_BASES;
        setRefBases(refBases);

        // test 1: 4xA -> 3xA contraction on left, 1xG ins on right, ie a left shift
        SimpleVariant variant = new SimpleVariant(CHR_1, 23, varRefBases, varAltBases);

        UltimaQualModel model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertEquals(BASE_SHIFT, model.type());

        UltimaBaseShift baseShiftModel = (UltimaBaseShift)model;
        assertNotNull(baseShiftModel.leftAdjust());
        assertNotNull(baseShiftModel.rightAdjust());
        assertTrue(baseShiftModel.isLeftShift());

        String readBases = BUFFER_REF_BASES + leftRefBases + varAltBases + rightRefBases;
        byte[] baseQualities = buildDefaultBaseQuals(readBases.length());
        byte[] t0Values = buildDefaultBaseQuals(readBases.length());
        byte[] tpValues = new byte[readBases.length()];

        SAMRecord read = buildUltimaRead(readBases, 1, baseQualities, tpValues, t0Values);

        byte calcQual = model.calculateQual(read, 22);
        assertEquals(50, calcQual);

        // test 2: 1xG ins on left, 4xT -> 3xT contraction on right

        // ref bases:   CAAA CCT TTTC
        // read bases:  CAAA GCC TTTC
        leftRefBases = "CAAA";
        varRefBases = "CCT";
        varAltBases = "GCC";
        rightRefBases = "TTTC";
        refBases = BUFFER_REF_BASES + leftRefBases + varRefBases + rightRefBases + BUFFER_REF_BASES;
        setRefBases(refBases);

        variant = new SimpleVariant(CHR_1, 23, varRefBases, varAltBases);

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertEquals(BASE_SHIFT, model.type());

        baseShiftModel = (UltimaBaseShift)model;
        assertNotNull(baseShiftModel.leftAdjust());
        assertNotNull(baseShiftModel.rightAdjust());
        assertTrue(baseShiftModel.isRightShift());

        // test 3
        // ref: ACGT ACC TTTC
        // alt: ACGT CCG TTTC - 1xA del on left, 1xG ins on right, left shift

        leftRefBases = "ACGT";
        varRefBases = "ACC";
        varAltBases = "CCG";
        rightRefBases = "TTTC";
        refBases = BUFFER_REF_BASES + leftRefBases + varRefBases + rightRefBases + BUFFER_REF_BASES;
        setRefBases(refBases);

        variant = new SimpleVariant(CHR_1, 23, varRefBases, varAltBases);

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertEquals(BASE_SHIFT, model.type());

        baseShiftModel = (UltimaBaseShift)model;
        assertNotNull(baseShiftModel.leftDeletion());
        assertNotNull(baseShiftModel.rightAdjust());
        assertTrue(baseShiftModel.isLeftShift());

        // test 4
        // ref: TAGG AT GGAA
        // alt: TAGG TG GGAA - 1xA del on left, 2xG -> 3xG on right, left shift
        leftRefBases = "TAGG";
        varRefBases = "AT";
        varAltBases = "TG";
        rightRefBases = "GGAA";
        refBases = BUFFER_REF_BASES + leftRefBases + varRefBases + rightRefBases + BUFFER_REF_BASES;
        setRefBases(refBases);

        variant = new SimpleVariant(CHR_1, 23, varRefBases, varAltBases);

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertEquals(BASE_SHIFT, model.type());

        baseShiftModel = (UltimaBaseShift)model;
        assertNotNull(baseShiftModel.leftDeletion());
        assertNotNull(baseShiftModel.rightAdjust());
        assertTrue(baseShiftModel.isLeftShift());

        // test 5:
        // ref: TAGG AT GGAA
        // alt: TAGG CA GGAA - 1xC ins on left, 1xT del on right
        varRefBases = "AT";
        varAltBases = "CA";
        refBases = BUFFER_REF_BASES + leftRefBases + varRefBases + rightRefBases + BUFFER_REF_BASES;
        setRefBases(refBases);

        variant = new SimpleVariant(CHR_1, 23, varRefBases, varAltBases);

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertEquals(BASE_SHIFT, model.type());

        baseShiftModel = (UltimaBaseShift)model;
        assertNotNull(baseShiftModel.leftAdjust());
        assertNotNull(baseShiftModel.rightDeletion());
        assertTrue(baseShiftModel.isRightShift());

        // test 6:
        // ref: TAGG AT GGAA
        // alt: TAGG GA GGAA - 2xG -> 3xG ins on left, 1xT del on right
        varRefBases = "AT";
        varAltBases = "GA";
        refBases = BUFFER_REF_BASES + leftRefBases + varRefBases + rightRefBases + BUFFER_REF_BASES;
        setRefBases(refBases);

        variant = new SimpleVariant(CHR_1, 23, varRefBases, varAltBases);

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertEquals(BASE_SHIFT, model.type());

        baseShiftModel = (UltimaBaseShift)model;
        assertNotNull(baseShiftModel.leftAdjust());
        assertNotNull(baseShiftModel.rightDeletion());
        assertTrue(baseShiftModel.isRightShift());

        // test 7: invalid
        varRefBases = "AT";
        varAltBases = "CG";
        refBases = BUFFER_REF_BASES + leftRefBases + varRefBases + rightRefBases + BUFFER_REF_BASES;
        setRefBases(refBases);

        variant = new SimpleVariant(CHR_1, 23, varRefBases, varAltBases);
        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertEquals(OTHER, model.type());

        // test 8: invalid
        varRefBases = "ATC";
        varAltBases = "TAG";
        refBases = BUFFER_REF_BASES + leftRefBases + varRefBases + rightRefBases + BUFFER_REF_BASES;
        setRefBases(refBases);

        variant = new SimpleVariant(CHR_1, 23, varRefBases, varAltBases);

        model = mModelBuilder.buildContext(variant, getStraddlingBases(refBases, variant));
        assertEquals(OTHER, model.type());
    }

    private static SAMRecord buildUltimaRead(
            final String readBases, final int readStart, final byte[] qualities, final byte[] tpValues, final byte[] t0Values)
    {
        return SeqTechTestUtils.buildUltimaRead(READ_ID_GENERATOR.nextId(), CHR_1, readStart, readBases, qualities, tpValues, t0Values);
    }

    private static byte[] getStraddlingBases(final String refBases, final SimpleVariant variant)
    {
        int refVarIndex = variant.position() - 1;

        byte[] straddlingBases = new byte[2];

        if(variant.isInsert())
            return straddlingBases;

        straddlingBases[0] = (byte)refBases.charAt(refVarIndex - 1);

        if(variant.isDelete())
        {
            straddlingBases[1] = (byte)refBases.charAt(refVarIndex + variant.refLength());
        }
        else
        {
            straddlingBases[1] = (byte)refBases.charAt(refVarIndex + variant.altLength());
        }

        return straddlingBases;
    }
}
