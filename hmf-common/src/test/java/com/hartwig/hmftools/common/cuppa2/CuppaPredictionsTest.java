package com.hartwig.hmftools.common.cuppa2;
import com.google.common.io.Resources;
import org.junit.Test;
import java.io.IOException;

public class CuppaPredictionsTest
{
    public CuppaPredictionsTest() throws IOException {}

    private static final String CUPPA_VIS_DATA_TSV = Resources.getResource("cuppa/cuppa_vis_data.tsv").getPath();
    CuppaPredictions cuppaPredictions = CuppaPredictions.fromTsv(CUPPA_VIS_DATA_TSV);

    @Test
    public void testPrintPredictions()
    {
        // cuppaPredictions.printPredictions(20);
    }

    @Test
    public void subsetByDataTypeReturnsCorrectDataType()
    {
        CuppaPredictionEntry predictionEntry = cuppaPredictions.subsetByDataType(Categories.DataType.FEAT_CONTRIB).get(0);
        assert predictionEntry.DataType.equals(Categories.DataType.FEAT_CONTRIB);
    }


    @Test
    public void sortByRankHasCorrectOrder()
    {
        CuppaPredictions sortedCuppaPredictions = cuppaPredictions.sortByRank();

        CuppaPredictionEntry entry1 = sortedCuppaPredictions.get(0);
        CuppaPredictionEntry entry2 = sortedCuppaPredictions.get(1);
        CuppaPredictionEntry entry3 = sortedCuppaPredictions.get(2);

        assert entry1.RankGroup == 0;
        assert entry1.Rank == 1;

        assert entry2.RankGroup == 0;
        assert entry2.Rank == 2;

        assert entry3.RankGroup == 0;
        assert entry3.Rank == 3;
    }
}

