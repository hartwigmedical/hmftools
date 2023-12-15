package com.hartwig.hmftools.orange.algo.cuppa2;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictions;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableCuppaPredictions;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableFeatureContributionEntry;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableProbabilityEntry;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableSignatureQuantileEntry;

import org.junit.Test;

public class CuppaPredictionsFactoryTest
{
    public CuppaPredictionsFactoryTest() throws IOException {}

    private static final String CUPPA_VIS_DATA_TSV = Resources.getResource("cuppa/cuppa_vis_data.tsv").getPath();
    CuppaPredictions cuppaPredictions = CuppaPredictions.fromTsv(CUPPA_VIS_DATA_TSV);

    @Test
    public void topPredictionIsCorrect()
    {
        ImmutableProbabilityEntry probabilityEntry = CuppaPredictionsFactory.getTopPrediction(cuppaPredictions);
        assert probabilityEntry.cancerType().equals("Breast: Triple negative");
    }

    @Test
    public void getProbabilitiesSucceeds()
    {
        ImmutableProbabilityEntry probabilityEntry = CuppaPredictionsFactory.getProbabilities(cuppaPredictions).get(0);
        // System.out.println(probabilityEntry.dataValue());
        assert true;
    }

    @Test
    public void getFeatureContributionEntriesSucceeds()
    {
        ImmutableFeatureContributionEntry featureContributionEntry = CuppaPredictionsFactory.getFeatureContributions(cuppaPredictions).get(0);
        // System.out.println(featureContributionEntry.dataValue());
        assert true;
    }

    @Test
    public void getSignatureQuantilesSucceeds()
    {
        ImmutableSignatureQuantileEntry signatureQuantileEntry = CuppaPredictionsFactory.getSignatureQuantiles(cuppaPredictions).get(0);
        // System.out.println(signatureQuantileEntry.dataValue());
        assert true;
    }

    @Test
    public void createImmutableDataSucceeds()
    {
        ImmutableCuppaPredictions data = CuppaPredictionsFactory.create(cuppaPredictions);
        assert data.topPrediction().cancerType().equals("Breast: Triple negative");
        assert data.probs().size()==320;
        assert data.featContribs().size()==640;
        assert data.sigQuantiles().size()==200;
    }
}
