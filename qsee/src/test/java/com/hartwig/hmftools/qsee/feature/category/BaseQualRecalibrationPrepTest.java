package com.hartwig.hmftools.qsee.feature.category;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;

import static com.hartwig.hmftools.qsee.feature.category.BaseQualRecalibrationPrep.aggregateChangeInQualPerOriginalQual;
import static com.hartwig.hmftools.qsee.feature.category.BaseQualRecalibrationPrep.aggregateChangeInQualPerTrinucContext;

import java.util.List;

import com.hartwig.hmftools.common.bam.ConsensusType;

import com.hartwig.hmftools.common.redux.BqrKey;
import com.hartwig.hmftools.common.redux.BqrRecord;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.category.bqr.BaseQualBin;
import com.hartwig.hmftools.qsee.feature.category.bqr.BaseQualBinner;
import com.hartwig.hmftools.qsee.feature.category.bqr.StandardBqrRecord;

import org.junit.Test;

public class BaseQualRecalibrationPrepTest
{
    private static final List<StandardBqrRecord> BQR_RECORDS = List.of(
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
            createBqrRecord(ConsensusType.NONE, 'G', 'T', "TGT", 10, (byte) 42, 41.5),

            // T>A_ATA: count below min count
            createBqrRecord(ConsensusType.NONE, 'T', 'A', "ATA", 2, (byte) 37, 39.0)
    );

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
    public void canAggregateChangeInQualPerTrinucContext()
    {
        BaseQualBinner baseQualBinner = new BaseQualBinner(SequencingType.ILLUMINA);
        List<Feature> features = aggregateChangeInQualPerTrinucContext(BQR_RECORDS, baseQualBinner);

        assertEquals(5, features.size());

        Feature actualFeature;

        actualFeature = features.get(0);
        assertEquals("ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:CCC;ReadOrientation:FORWARD;OriginalQualBin:HIGH (30+)", actualFeature.key().name());
        assertEquals(2.0, actualFeature.value(), 0.001);

        actualFeature = features.get(1);
        assertEquals("ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:CCC;ReadOrientation:REVERSE;OriginalQualBin:HIGH (30+)", actualFeature.key().name());
        assertEquals(-0.5, actualFeature.value(), 0.001);

        actualFeature = features.get(2);
        assertEquals("ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:ACA;ReadOrientation:FORWARD;OriginalQualBin:HIGH (30+)", actualFeature.key().name());
        assertEquals(4.0, actualFeature.value(), 0.001);

        actualFeature = features.get(3);
        assertEquals("ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:ACA;ReadOrientation:REVERSE;OriginalQualBin:HIGH (30+)", actualFeature.key().name());
        assertEquals(-0.5, actualFeature.value(), 0.001);

        actualFeature = features.get(4);
        assertEquals("ReadType:NONE;StandardMutation:T>A;StandardTrinucContext:ATA;ReadOrientation:FORWARD;OriginalQualBin:HIGH (30+)", actualFeature.key().name());
        assertTrue(Double.isNaN(actualFeature.value()));
    }

    @Test
    public void canAggregateChangeInQualPerOriginalQual()
    {
        BaseQualBinner baseQualBinner = new BaseQualBinner(SequencingType.ILLUMINA);
        List<Feature> features = aggregateChangeInQualPerOriginalQual(BQR_RECORDS, baseQualBinner);

        assertEquals(5, features.size());

        Feature actualFeature;

        actualFeature = features.get(0);
        assertEquals("ReadType:DUAL;StandardMutation:C>G;StandardTrinucContext:ANY;ReadOrientation:FORWARD;OriginalQualBin:LOW (0-29)", actualFeature.key().name());
        assertEquals(2.0, actualFeature.value(), 0.001);

        actualFeature = features.get(1);
        assertEquals("ReadType:DUAL;StandardMutation:C>G;StandardTrinucContext:ANY;ReadOrientation:REVERSE;OriginalQualBin:LOW (0-29)", actualFeature.key().name());
        assertEquals(-0.5, actualFeature.value(), 0.001);

        actualFeature = features.get(2);
        assertEquals("ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:ANY;ReadOrientation:FORWARD;OriginalQualBin:HIGH (30+)", actualFeature.key().name());
        assertEquals(3.118, actualFeature.value(), 0.001);

        actualFeature = features.get(3);
        assertEquals("ReadType:NONE;StandardMutation:C>A;StandardTrinucContext:ANY;ReadOrientation:REVERSE;OriginalQualBin:HIGH (30+)", actualFeature.key().name());
        assertEquals(-0.5, actualFeature.value(), 0.001);

        actualFeature = features.get(4);
        assertEquals("ReadType:NONE;StandardMutation:T>A;StandardTrinucContext:ANY;ReadOrientation:FORWARD;OriginalQualBin:HIGH (30+)", actualFeature.key().name());
        assertTrue(Double.isNaN(actualFeature.value()));
    }

    private static StandardBqrRecord createBqrRecord(ConsensusType readType, char refBase, char altBase, String trinucContext,
            int count, int originalQual, double recalibratedQual)
    {
        BqrKey key = new BqrKey((byte) refBase, (byte) altBase, trinucContext.getBytes(), (byte) originalQual, readType);
        return StandardBqrRecord.from(new BqrRecord(key, count, recalibratedQual));
    }
}
