package com.hartwig.hmftools.orange.algo.cuppa2;
import java.io.IOException;

import com.hartwig.hmftools.common.cuppa2.CuppaPredictions;
import com.hartwig.hmftools.common.cuppa2.FeatureContributionEntry;
import com.hartwig.hmftools.common.cuppa2.ProbabilityEntry;
import com.hartwig.hmftools.common.cuppa2.SignatureQuantileEntry;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableCuppaPredictionData;

import org.junit.Test;

public class CuppaPredictionDataFactoryTest
{
    public CuppaPredictionDataFactoryTest() throws IOException {}

    private static final String CUPPA_VIS_DATA_TSV = "/Users/lnguyen/Hartwig/hartwigmedical/hmftools/hmf-common/src/test/resources/cuppa/cuppa_vis_data.tsv";
    CuppaPredictions cuppaPredictions = CuppaPredictions.fromTsv(CUPPA_VIS_DATA_TSV);

    @Test
    public void topPredictionIsCorrect()
    {
        ProbabilityEntry probabilityEntry = CuppaPredictionDataFactory.getTopPrediction(cuppaPredictions);
        assert probabilityEntry.CancerType.equals("Breast: Triple negative");
    }

    @Test
    public void getProbabilityEntriesSucceeds()
    {
        ProbabilityEntry probabilityEntry = CuppaPredictionDataFactory.getProbabilities(cuppaPredictions).get(0);
        // System.out.println(probabilityEntry.DataValue);
        assert true;
    }

    @Test
    public void getFeatureContributionEntriesSucceeds()
    {
        FeatureContributionEntry featureContributionEntry = CuppaPredictionDataFactory.getFeatureContributions(cuppaPredictions).get(0);
        // System.out.println(featureContributionEntry.DataValue);
        assert true;
    }

    @Test
    public void getSignatureQuantilesSucceeds()
    {
        SignatureQuantileEntry signatureQuantileEntry = CuppaPredictionDataFactory.getSignatureQuantiles(cuppaPredictions).get(0);
        // System.out.println(signatureQuantileEntry.DataValue);
        assert true;
    }

    @Test
    public void createImmutableDataSucceeds()
    {
        ImmutableCuppaPredictionData data = CuppaPredictionDataFactory.create(cuppaPredictions);
        assert data.topPrediction().CancerType.equals("Breast: Triple negative");
        assert data.probs().size()==320;
        assert data.featContribs().size()==640;
        assert data.sigQuantiles().size()==200;
    }
}
