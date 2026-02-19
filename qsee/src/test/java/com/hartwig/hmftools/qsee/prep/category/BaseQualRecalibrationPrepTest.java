package com.hartwig.hmftools.qsee.prep.category;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import static com.hartwig.hmftools.qsee.prep.category.BaseQualRecalibrationPrep.calcChangeInQualPerOriginalQual;
import static com.hartwig.hmftools.qsee.prep.category.BaseQualRecalibrationPrep.calcChangeInQualPerTrinucContext;

import java.util.List;

import com.hartwig.hmftools.common.bam.ConsensusType;

import com.hartwig.hmftools.common.redux.BqrKey;
import com.hartwig.hmftools.common.redux.BqrRecord;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.prep.category.bqr.BaseQualBin;
import com.hartwig.hmftools.qsee.prep.category.bqr.BaseQualBinner;

import org.junit.Test;

public class BaseQualRecalibrationPrepTest
{
    private static final List<BqrRecord> BQR_RECORDS = List.of(
            // C>G_CCC
            createBqrRecord(ConsensusType.DUAL, 'C', 'G', "CCC", 90, (byte) 11, 13.0),
            createBqrRecord(ConsensusType.DUAL, 'C', 'G', "CCC", 90, (byte) 17, 19.0),
            createBqrRecord(ConsensusType.DUAL, 'G', 'C', "GGG", 10, (byte) 11, 10.5),
            createBqrRecord(ConsensusType.DUAL, 'G', 'C', "GGG", 10, (byte) 17, 16.5),

            // C>A_CCC
            createBqrRecord(ConsensusType.NONE, 'C', 'A', "CCC", 70, (byte) 37, 39.0),
            createBqrRecord(ConsensusType.NONE, 'C', 'A', "CCC", 70, (byte) 42, 44.0),
            createBqrRecord(ConsensusType.NONE, 'G', 'T', "GGG", 30, (byte) 37, 36.5),
            createBqrRecord(ConsensusType.NONE, 'G', 'T', "GGG", 30, (byte) 42, 41.5),

            // C>A_ACA
            createBqrRecord(ConsensusType.NONE, 'C', 'A', "ACA", 90, (byte) 37, 41.0),
            createBqrRecord(ConsensusType.NONE, 'C', 'A', "ACA", 90, (byte) 42, 46.0),
            createBqrRecord(ConsensusType.NONE, 'G', 'T', "TGT", 10, (byte) 37, 36.5),
            createBqrRecord(ConsensusType.NONE, 'G', 'T', "TGT", 10, (byte) 42, 41.5)
    );

    @Test
    public void canStandardiseBases()
    {
        // Already standard: should remain unchanged
        BqrRecord record1 = createBqrRecord(ConsensusType.NONE, 'C', 'A', "ACC",
                0, 0, 0);

        assertEquals('C', record1.Key.Ref);
        assertEquals('A', record1.Key.Alt);
        assertEquals("ACC", new String(record1.Key.TrinucleotideContext));

        // Not standard: should be reverse complemented
        BqrRecord record2 = createBqrRecord(ConsensusType.NONE, 'G', 'T', "GGA",
                0, 0, 0);

        assertEquals('C', record2.Key.Ref);
        assertEquals('A', record2.Key.Alt);
        assertEquals("TCC", new String(record2.Key.TrinucleotideContext));
    }

    @Test
    public void canGetBaseQualBinRangesForIllumina()
    {
        BaseQualBinner baseQualBinner = new BaseQualBinner(SequencingType.ILLUMINA);

        assertEquals(0, baseQualBinner.binRanges().get(BaseQualBin.LOW).lowerBound());
        assertEquals(29, baseQualBinner.binRanges().get(BaseQualBin.LOW).upperBound());

        assertFalse(baseQualBinner.binRanges().containsKey(BaseQualBin.MEDIUM));

        assertEquals(30, baseQualBinner.binRanges().get(BaseQualBin.HIGH).lowerBound());
        assertEquals(Byte.MAX_VALUE, baseQualBinner.binRanges().get(BaseQualBin.HIGH).upperBound());
    }

    @Test
    public void canCalcChangeInQualPerTrinucContext()
    {
        BaseQualBinner baseQualBinner = new BaseQualBinner(SequencingType.ILLUMINA);
        List<Feature> features = calcChangeInQualPerTrinucContext(BQR_RECORDS, baseQualBinner);

        assertEquals(2, features.size());

        Feature actualFeature;

        actualFeature = features.get(0);
        assertEquals("ReadType=NONE;StandardMutation=C>A;StandardTrinucContext=CCC;OriginalQualBin=HIGH (30+)", actualFeature.key().name());
        assertEquals(1.205, actualFeature.value(), 0.001);

        actualFeature = features.get(1);
        assertEquals("ReadType=NONE;StandardMutation=C>A;StandardTrinucContext=ACA;OriginalQualBin=HIGH (30+)", actualFeature.key().name());
        assertEquals(3.386, actualFeature.value(), 0.001);
    }

    @Test
    public void canCalcChangeInQualPerOriginalQual()
    {
        BaseQualBinner baseQualBinner = new BaseQualBinner(SequencingType.ILLUMINA);
        List<Feature> features = calcChangeInQualPerOriginalQual(BQR_RECORDS, baseQualBinner);

        assertEquals(2, features.size());

        Feature actualFeature;

        actualFeature = features.get(0);
        assertEquals("ReadType=DUAL;StandardMutation=C>G;OriginalQualBin=LOW (0-29)", actualFeature.key().name());
        assertEquals(1.659, actualFeature.value(), 0.001);

        actualFeature = features.get(1);
        assertEquals("ReadType=NONE;StandardMutation=C>A;OriginalQualBin=HIGH (30+)", actualFeature.key().name());
        assertEquals(2.295, actualFeature.value(), 0.001);
    }

    private static BqrRecord createBqrRecord(ConsensusType readType, char refBase, char altBase, String trinucContext,
            int count, int originalQual, double recalibratedQual)
    {
        BqrKey key = new BqrKey((byte) refBase, (byte) altBase, trinucContext.getBytes(), (byte) originalQual, readType);
        BqrRecord record = new BqrRecord(key, count, recalibratedQual);
        record = BaseQualRecalibrationPrep.standardiseBases(record);

        return record;
    }
}
