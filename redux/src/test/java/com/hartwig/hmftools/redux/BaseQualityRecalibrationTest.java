package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.redux.bqr.BaseQualRecalibration.convertToRecords;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.redux.BqrKey;
import com.hartwig.hmftools.common.redux.BqrRecord;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.redux.bqr.BaseQualityData;
import com.hartwig.hmftools.redux.bqr.BaseQualityResults;
import com.hartwig.hmftools.redux.bqr.BqrKeyCounter;
import com.hartwig.hmftools.redux.bqr.BqrRegionReader;

import org.junit.Test;

public class BaseQualityRecalibrationTest
{
    private static final MockRefGenome REF_GENOME = new MockRefGenome(false);

    static
    {
        REF_GENOME.RefGenomeMap.put(CHR_1, generateRandomBases(300));
    }

    @Test
    public void testBaseQualityCounts()
    {
        BqrRegionReader bqrCounter = new BqrRegionReader(null, new BaseQualityResults(), Collections.emptyList());

        bqrCounter.initialise(new ChrBaseRegion(CHR_1, 100, 250), REF_GENOME);

        int pos1 = 100;
        BqrKey key1 = createKey('A', 'G', 30, ConsensusType.NONE);
        BqrKey key2 = createKey('A', 'A', 20, ConsensusType.NONE);
        BqrKey key3 = createKey('A', 'G', 15, ConsensusType.NONE); // a repeated alt

        BaseQualityData bqData1 = bqrCounter.getOrCreateBaseQualData(pos1, key1.Ref, key1.TrinucleotideContext);

        bqData1.processReadBase(ConsensusType.NONE, key1.Alt, key1.Quality, true);

        for(int i = 0; i < 10; ++i)
        {
            bqData1.processReadBase(ConsensusType.NONE, key2.Alt, key2.Quality, true);
        }

        for(int i = 0; i < 3; ++i)
        {
            bqData1.processReadBase(ConsensusType.NONE, key3.Alt, key3.Quality, true);
        }

        // repeated alt at different locations
        int pos2 = 150;
        BqrKey key4 = createKey('C', 'G', 25, ConsensusType.NONE); // another repeated alt
        BaseQualityData bqData2 = bqrCounter.getOrCreateBaseQualData(pos2, key4.Ref, key4.TrinucleotideContext);

        for(int i = 0; i < 4; ++i)
        {
            bqData2.processReadBase(ConsensusType.NONE, key4.Alt, key4.Quality, true);
        }

        int pos3 = 200;
        BqrKey key5 = createKey('A', 'G', 20, ConsensusType.NONE); // an alt but not repeated
        BaseQualityData bqData3 = bqrCounter.getOrCreateBaseQualData(pos3, key5.Ref, key5.TrinucleotideContext);

        bqData3.processReadBase(ConsensusType.NONE, key5.Alt, key5.Quality, true);

        for(int i = 0; i < 9; ++i)
        {
            bqData3.processReadBase(ConsensusType.NONE, key2.Ref, key2.Quality, true); // AF of 10% but count of 1 is permitted
        }

        bqrCounter.buildQualityCounts();

        Collection<BqrKeyCounter> qualityCounts = bqrCounter.getQualityCounts();

        BqrKeyCounter qc = qualityCounts.stream().filter(x -> x.Key.equals(key1)).findFirst().orElse(null);
        assertTrue(qc == null);

        qc = qualityCounts.stream().filter(x -> x.Key.equals(key2)).findFirst().orElse(null);
        assertTrue(qc != null);
        assertEquals(19, qc.count());

        qc = qualityCounts.stream().filter(x -> x.Key.equals(key3)).findFirst().orElse(null);
        assertTrue(qc == null);

        qc = qualityCounts.stream().filter(x -> x.Key.equals(key4)).findFirst().orElse(null);
        assertTrue(qc == null);

        qc = qualityCounts.stream().filter(x -> x.Key.equals(key5)).findFirst().orElse(null);
        assertTrue(qc != null);
        assertEquals(1, qc.count());
    }

    @Test
    public void testBaseQualityCountsByReadType()
    {
        BqrRegionReader bqrCounter = new BqrRegionReader(null, new BaseQualityResults(), Collections.emptyList());

        bqrCounter.initialise(new ChrBaseRegion(CHR_1, 100, 250), REF_GENOME);

        char refBase = 'G';
        char altBase = 'G';
        byte qual1 = 37;
        byte qual2 = 25;

        addReadBaseQual(bqrCounter, 100, refBase, altBase, ConsensusType.NONE, qual1);
        addReadBaseQual(bqrCounter, 100, refBase, altBase, ConsensusType.NONE, qual1);
        addReadBaseQual(bqrCounter, 100, refBase, altBase, ConsensusType.NONE, qual1);

        addReadBaseQual(bqrCounter, 100, refBase, altBase, ConsensusType.SINGLE, qual1);
        addReadBaseQual(bqrCounter, 100, refBase, altBase, ConsensusType.SINGLE, qual1);

        addReadBaseQual(bqrCounter, 100, refBase, altBase, ConsensusType.DUAL, qual1);
        addReadBaseQual(bqrCounter, 101, refBase, altBase, ConsensusType.DUAL, qual1);
        addReadBaseQual(bqrCounter, 102, refBase, altBase, ConsensusType.DUAL, qual1);
        addReadBaseQual(bqrCounter, 103, refBase, altBase, ConsensusType.DUAL, qual1);

        addReadBaseQual(bqrCounter, 101, refBase, altBase, ConsensusType.DUAL, qual2);
        addReadBaseQual(bqrCounter, 102, refBase, altBase, ConsensusType.DUAL, qual2);

        BqrKey keyNone = createKey(refBase, altBase, qual1, ConsensusType.NONE);
        BqrKey keySingle = createKey(refBase, altBase, qual1, ConsensusType.SINGLE);
        BqrKey keyDualQ1 = createKey(refBase, altBase, qual1, ConsensusType.DUAL);
        BqrKey keyDualQ2 = createKey(refBase, altBase, qual2, ConsensusType.DUAL);

        bqrCounter.buildQualityCounts();

        Collection<BqrKeyCounter> qualityCounts = bqrCounter.getQualityCounts();

        BqrKeyCounter qc = qualityCounts.stream().filter(x -> x.Key.equals(keyNone)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(3, qc.count());

        qc = qualityCounts.stream().filter(x -> x.Key.equals(keySingle)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(2, qc.count());

        qc = qualityCounts.stream().filter(x -> x.Key.equals(keyDualQ1)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(4, qc.count());

        qc = qualityCounts.stream().filter(x -> x.Key.equals(keyDualQ2)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(2, qc.count());
    }

    private static void addReadBaseQual(
            final BqrRegionReader bqrCounter, int position, char ref, char alt, final ConsensusType readType, byte quality)
    {
        byte[] context = new byte[] { DNA_BASE_BYTES[0], (byte)ref, DNA_BASE_BYTES[0]};
        BaseQualityData baseQualityData = bqrCounter.getOrCreateBaseQualData(position, (byte)ref, context);
        baseQualityData.processReadBase(readType, (byte)alt, quality, true);
    }

    private BqrKey createKey(char ref, char alt, int qual, final ConsensusType readType)
    {
        byte[] context = new byte[] { 65,  (byte)ref, 65};
        return new BqrKey((byte)ref, (byte)alt, context, (byte)qual, readType);
    }

    @Test
    public void testBaseQualityAltVafChecks()
    {
        BqrRegionReader bqrCounter = new BqrRegionReader(null, new BaseQualityResults(), Collections.emptyList());

        bqrCounter.initialise(new ChrBaseRegion(CHR_1, 1, 10), REF_GENOME);

        byte baseQual = 37;
        int position = 5;

        byte ref = DNA_BASE_BYTES[1];
        byte alt = DNA_BASE_BYTES[2];
        byte[] tnContext = new byte[] { DNA_BASE_BYTES[0], ref, DNA_BASE_BYTES[0]};
        BaseQualityData baseQualityData = bqrCounter.getOrCreateBaseQualData(position, ref, tnContext);

        // both read types pass the VAF / AD tests

        boolean posStrand = true;

        addReadCounts(baseQualityData, ConsensusType.NONE, ref, baseQual, posStrand, 100);
        addReadCounts(baseQualityData, ConsensusType.DUAL, ref, baseQual, posStrand, 100);

        addReadCounts(baseQualityData, ConsensusType.NONE, alt, baseQual, posStrand, 2);
        addReadCounts(baseQualityData, ConsensusType.DUAL, alt, baseQual, posStrand, 2);

        bqrCounter.buildQualityCounts();

        Collection<BqrKeyCounter> qualityCounts = bqrCounter.getQualityCounts();

        BqrKey keyNoneRef = new BqrKey(ref, ref, tnContext, baseQual, ConsensusType.NONE);
        BqrKey keyNoneAlt = new BqrKey(ref, alt, tnContext, baseQual, ConsensusType.NONE);
        BqrKey keyDualRef = new BqrKey(ref, ref, tnContext, baseQual, ConsensusType.DUAL);
        BqrKey keyDualAlt = new BqrKey(ref, alt, tnContext, baseQual, ConsensusType.DUAL);

        BqrKeyCounter qc = qualityCounts.stream().filter(x -> x.Key.equals(keyNoneRef)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(100, qc.count());

        qc = qualityCounts.stream().filter(x -> x.Key.equals(keyDualRef)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(100, qc.count());

        qc = qualityCounts.stream().filter(x -> x.Key.equals(keyNoneAlt)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(2, qc.count());

        qc = qualityCounts.stream().filter(x -> x.Key.equals(keyDualAlt)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(2, qc.count());

        // the standard read types will exceed the VAF tests but the dual will not, and so both are dropped
        bqrCounter.initialise(new ChrBaseRegion(CHR_1, 1, 10), REF_GENOME);

        baseQualityData = bqrCounter.getOrCreateBaseQualData(position, ref, tnContext);

        addReadCounts(baseQualityData, ConsensusType.NONE, ref, baseQual, posStrand, 100);
        addReadCounts(baseQualityData, ConsensusType.DUAL, ref, baseQual, posStrand, 100);

        addReadCounts(baseQualityData, ConsensusType.NONE, alt, baseQual, posStrand, 6);
        addReadCounts(baseQualityData, ConsensusType.DUAL, alt, baseQual, posStrand, 2);

        bqrCounter.buildQualityCounts();

        qualityCounts = bqrCounter.getQualityCounts();

        qc = qualityCounts.stream().filter(x -> x.Key.equals(keyNoneAlt)).findFirst().orElse(null);
        assertNull(qc);
        assertNull(qc);

        // the dual will fail and the standard will pass
        bqrCounter.initialise(new ChrBaseRegion(CHR_1, 1, 10), REF_GENOME);

        baseQualityData = bqrCounter.getOrCreateBaseQualData(position, ref, tnContext);

        addReadCounts(baseQualityData, ConsensusType.NONE, ref, baseQual, posStrand, 100);
        addReadCounts(baseQualityData, ConsensusType.DUAL, ref, baseQual, posStrand, 100);

        addReadCounts(baseQualityData, ConsensusType.NONE, alt, baseQual, posStrand, 2);
        addReadCounts(baseQualityData, ConsensusType.DUAL, alt, baseQual, posStrand, 3);

        bqrCounter.buildQualityCounts();

        qualityCounts = bqrCounter.getQualityCounts();

        qc = qualityCounts.stream().filter(x -> x.Key.equals(keyNoneAlt)).findFirst().orElse(null);
        assertNull(qc);
        assertNull(qc);
    }

    @Test
    public void testBaseQualityReadStrands()
    {
        BqrRegionReader bqrCounter = new BqrRegionReader(null, new BaseQualityResults(), Collections.emptyList());

        bqrCounter.initialise(new ChrBaseRegion(CHR_1, 1, 10), REF_GENOME);

        byte baseQual = 37;
        int position = 3;

        byte ref = DNA_BASE_BYTES[1];
        byte alt = DNA_BASE_BYTES[2];
        byte[] tnContext = new byte[] { DNA_BASE_BYTES[0], ref, DNA_BASE_BYTES[0] };
        BaseQualityData baseQualityData = bqrCounter.getOrCreateBaseQualData(position, ref, tnContext);

        addReadCounts(baseQualityData, ConsensusType.NONE, ref, baseQual, true, 50);
        addReadCounts(baseQualityData, ConsensusType.NONE, ref, baseQual, false, 50);

        addReadCounts(baseQualityData, ConsensusType.NONE, alt, baseQual, true, 1);
        addReadCounts(baseQualityData, ConsensusType.NONE, alt, baseQual, false, 1);

        byte[] tnContextReversed = Nucleotides.reverseComplementBases(tnContext);
        byte refReversed = Nucleotides.swapDnaBase(ref);
        byte altReversed = Nucleotides.swapDnaBase(alt);

        BaseQualityData baseQualityDataRev = bqrCounter.getOrCreateBaseQualData(position + 4, refReversed, tnContextReversed);

        addReadCounts(baseQualityDataRev, ConsensusType.NONE, refReversed, baseQual, true, 50);
        addReadCounts(baseQualityDataRev, ConsensusType.NONE, refReversed, baseQual, false, 50);
        addReadCounts(baseQualityDataRev, ConsensusType.NONE, altReversed, baseQual, true, 1);
        addReadCounts(baseQualityDataRev, ConsensusType.NONE, altReversed, baseQual, false, 1);

        bqrCounter.buildQualityCounts();

        Collection<BqrKeyCounter> qualityCounts = bqrCounter.getQualityCounts();

        BqrKey keyRef = new BqrKey(ref, ref, tnContext, baseQual, ConsensusType.NONE);
        BqrKey keyRefReversed = new BqrKey(refReversed, refReversed, tnContextReversed, baseQual, ConsensusType.NONE);

        BqrKey keyAlt = new BqrKey(ref, alt, tnContext, baseQual, ConsensusType.NONE);
        BqrKey keyAltReversed = new BqrKey(refReversed, altReversed, tnContextReversed, baseQual, ConsensusType.NONE);

        BqrKeyCounter qc = qualityCounts.stream().filter(x -> x.Key.equals(keyRef)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(100, qc.count());

        qc = qualityCounts.stream().filter(x -> x.Key.equals(keyRefReversed)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(100, qc.count());

        qc = qualityCounts.stream().filter(x -> x.Key.equals(keyAlt)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(2, qc.count());

        qc = qualityCounts.stream().filter(x -> x.Key.equals(keyAltReversed)).findFirst().orElse(null);
        assertNotNull(qc);
        assertEquals(2, qc.count());
    }

    private static void addReadCounts(
            final BaseQualityData baseQualityData, final ConsensusType readType, byte alt, byte baseQual, boolean posStrand, int count)
    {
        for(int i = 0; i < count; ++i)
        {
            baseQualityData.processReadBase(readType, alt, baseQual, posStrand);
        }
    }

    @Test
    public void testBaseQualRecalibrationCalcs()
    {
        Map<BqrKey,Integer> allQualityCounts = Maps.newHashMap();

        byte aBase = DNA_BASE_BYTES[0];
        byte cBase = DNA_BASE_BYTES[1];
        byte gBase = DNA_BASE_BYTES[2];

        byte qualHigh = 37;

        byte[] triNucContext1 = new byte[] {gBase, aBase, gBase};
        byte[] triNucContext2 = new byte[] {gBase, cBase, gBase};

        BqrKey aRefKey = new BqrKey(aBase, aBase, triNucContext1, qualHigh, ConsensusType.NONE);
        BqrKey cRefKey = new BqrKey(cBase, cBase, triNucContext2, qualHigh, ConsensusType.NONE);
        BqrKey cAltKey = new BqrKey(aBase, cBase, triNucContext1, qualHigh, ConsensusType.NONE);

        allQualityCounts.put(aRefKey, 1000);
        allQualityCounts.put(cRefKey, 4000);
        allQualityCounts.put(cAltKey, 10);

        byte qualLow = 10;
        BqrKey aRefKeyLow = new BqrKey(aBase, aBase, triNucContext1, qualLow, ConsensusType.NONE);
        BqrKey cRefKeyLow = new BqrKey(cBase, cBase, triNucContext2, qualLow, ConsensusType.NONE);
        BqrKey cAltKeyLow = new BqrKey(aBase, cBase, triNucContext1, qualLow, ConsensusType.NONE);

        allQualityCounts.put(aRefKeyLow, 2000);
        allQualityCounts.put(cRefKeyLow, 2000);
        allQualityCounts.put(cAltKeyLow, 2000);

        List<BqrRecord> bqrRecords = convertToRecords(allQualityCounts);

        assertEquals(16, bqrRecords.size());

        BqrRecord rec1 = bqrRecords.stream().filter(x -> x.Key.equals(cAltKey)).findFirst().orElse(null);
        assertEquals(25.2, rec1.RecalibratedQuality, 0.1);

        BqrKey aAltKey = new BqrKey(aBase, gBase, triNucContext1, qualHigh, ConsensusType.NONE);
        BqrRecord rec2 = bqrRecords.stream().filter(x -> x.Key.equals(aAltKey)).findFirst().orElse(null);
        assertEquals(37, rec2.RecalibratedQuality, 0.1);
    }
}
