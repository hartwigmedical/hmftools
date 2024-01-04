package com.hartwig.hmftools.common.cuppa2;
import com.google.common.io.Resources;
import org.junit.Test;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

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
    public void subsetByClfNameReturnsOneClfName()
    {
        List<CuppaPredictionEntry> predictionEntries = cuppaPredictions
                .subsetByClfName(Categories.ClfName.DNA_COMBINED)
                .PredictionEntries;

        List<Categories.ClfName> clfNames = predictionEntries.stream()
                .map(CuppaPredictionEntry::getClfName)
                .distinct()
                .collect(Collectors.toList());

        assert clfNames.size() == 1;
        assert clfNames.get(0).equals(Categories.ClfName.DNA_COMBINED);
    }

    @Test
    public void subsetByCancerTypeReturnsOneCancerType()
    {
        List<CuppaPredictionEntry> predictionEntries = cuppaPredictions
                .subsetByDataType(Categories.DataType.PROB)
                .subsetByCancerType("Breast: Triple negative")
                .PredictionEntries;

        List<String> cancerTypes = predictionEntries.stream()
                .map(CuppaPredictionEntry::getCancerType)
                .distinct().collect(Collectors.toList());

        assert cancerTypes.size() == 1;
        assert cancerTypes.get(0).equals("Breast: Triple negative");
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

