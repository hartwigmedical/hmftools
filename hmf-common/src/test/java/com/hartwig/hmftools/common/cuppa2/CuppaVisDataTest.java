package com.hartwig.hmftools.common.cuppa2;
import com.google.common.io.Resources;
import org.junit.Test;
import java.io.IOException;

public class CuppaVisDataTest
{
    public CuppaVisDataTest() throws IOException {}

    private static final String CUPPA_VIS_DATA_TSV = Resources.getResource("cuppa/cuppa_vis_data.tsv").getPath();
    CuppaVisData cuppaVisData = CuppaVisData.fromTsv(CUPPA_VIS_DATA_TSV);

    @Test
    public void subsetByDataTypeReturnsCorrectDataType()
    {
        CuppaVisDataEntry visDataEntry = cuppaVisData.subsetByDataType(Categories.DataType.FEAT_CONTRIB).get(0);
        assert visDataEntry.DataType.equals(Categories.DataType.FEAT_CONTRIB);
    }


    @Test
    public void sortByRankHasCorrectOrder()
    {
        CuppaVisData sortedCuppaVisData = cuppaVisData.sortByRank();

        CuppaVisDataEntry entry1 = sortedCuppaVisData.get(0);
        CuppaVisDataEntry entry2 = sortedCuppaVisData.get(1);
        CuppaVisDataEntry entry3 = sortedCuppaVisData.get(2);

        assert entry1.RankGroup == 0;
        assert entry1.Rank == 1;

        assert entry2.RankGroup == 0;
        assert entry2.Rank == 2;

        assert entry3.RankGroup == 0;
        assert entry3.Rank == 3;
    }
}

