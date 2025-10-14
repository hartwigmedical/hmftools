package prep;

import static org.junit.Assert.assertEquals;

import static prep.ReduxBqrPrep.calcChangeInQualPerOriginalQual;
import static prep.ReduxBqrPrep.calcChangeInQualPerTrinucContext;

import java.util.List;

import com.hartwig.hmftools.common.bam.ConsensusType;

import feature.FeatureValue;
import prep.ReduxBqrPrep.ExtendedBqrRecord;

import org.junit.Test;

public class ReduxBqrPrepTest
{
    private static final List<ExtendedBqrRecord> BQR_RECORDS = List.of(
            // C>G_CCC
            new ExtendedBqrRecord(ConsensusType.DUAL, 'C','G', "CCC", 90, (byte) 11, 13.0),
            new ExtendedBqrRecord(ConsensusType.DUAL, 'C','G', "CCC", 90, (byte) 17, 19.0),
            new ExtendedBqrRecord(ConsensusType.DUAL, 'G','C', "GGG", 10, (byte) 11, 10.5),
            new ExtendedBqrRecord(ConsensusType.DUAL, 'G','C', "GGG", 10, (byte) 17, 16.5),

            // C>A_CCC
            new ExtendedBqrRecord(ConsensusType.NONE, 'C','A', "CCC", 70, (byte) 37, 39.0),
            new ExtendedBqrRecord(ConsensusType.NONE, 'C','A', "CCC", 70, (byte) 42, 44.0),
            new ExtendedBqrRecord(ConsensusType.NONE, 'G','T', "GGG", 30, (byte) 37, 36.5),
            new ExtendedBqrRecord(ConsensusType.NONE, 'G','T', "GGG", 30, (byte) 42, 41.5),

            // C>A_ACA
            new ExtendedBqrRecord(ConsensusType.NONE, 'C','A', "ACA", 90, (byte) 37, 41.0),
            new ExtendedBqrRecord(ConsensusType.NONE, 'C','A', "ACA", 90, (byte) 42, 46.0),
            new ExtendedBqrRecord(ConsensusType.NONE, 'G','T', "TGT", 10, (byte) 37, 36.5),
            new ExtendedBqrRecord(ConsensusType.NONE, 'G','T', "TGT", 10, (byte) 42, 41.5)

        );

    @Test
    public void canCalcChangeInQualPerTrinucContext()
    {
        List<FeatureValue<Double>> featureValues = calcChangeInQualPerTrinucContext(BQR_RECORDS);

        assertEquals(2, featureValues.size());

        FeatureValue<Double> actualFeature;

        actualFeature = featureValues.get(0);
        assertEquals("readType=NONE;standardMutation=C>A;standardTrinucContext=CCC", actualFeature.mKey);
        assertEquals(1.25, actualFeature.mValue, 0.001);

        actualFeature = featureValues.get(1);
        assertEquals("readType=NONE;standardMutation=C>A;standardTrinucContext=ACA", actualFeature.mKey);
        assertEquals(3.55, actualFeature.mValue, 0.001);
    }

    @Test
    public void canCalcChangeInQualPerOriginalQual()
    {
        List<FeatureValue<Double>> featureValues = calcChangeInQualPerOriginalQual(BQR_RECORDS);

        assertEquals(2, featureValues.size());

        FeatureValue<Double> actualFeature;

        actualFeature = featureValues.get(0);
        assertEquals("readType=DUAL;standardMutation=C>G;originalQualBin=0-19", actualFeature.mKey);
        assertEquals(1.75, actualFeature.mValue, 0.001);

        actualFeature = featureValues.get(1);
        assertEquals("readType=NONE;standardMutation=C>A;originalQualBin=30+", actualFeature.mKey);
        assertEquals(2.40, actualFeature.mValue, 0.001);
    }
}
