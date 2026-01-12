package com.hartwig.hmftools.patientdb;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.cuppa.ClassifierGroup;
import com.hartwig.hmftools.common.cuppa.ClassifierName;
import com.hartwig.hmftools.common.cuppa.CuppaPredictionEntry;
import com.hartwig.hmftools.common.cuppa.CuppaPredictions;
import com.hartwig.hmftools.common.cuppa.DataType;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.CuppaRecord;

import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class LoadCuppaDataTest extends DatabaseAutoSetup
{
    private static final String TEST_SAMPLE_ID = "example";
    private static final int TOP_N_PROBS = 2;

    private static CuppaPredictions CUPPA_PREDICTIONS;
    private static List<CuppaRecord> CUPPA_RECORDS;

    @BeforeClass
    public static void loadInputData() throws Exception
    {
        List<CuppaPredictionEntry> predictionEntries = List.of(
                createTestPredictionEntry("CANCER_TYPE_1", 0.1, 2, 0),
                createTestPredictionEntry("CANCER_TYPE_2", 0.9, 1, 0),
                createTestPredictionEntry("CANCER_TYPE_3", 0.0, 3, 0),

                createTestPredictionEntry("CANCER_TYPE_1", 0.2, 2, 0),
                createTestPredictionEntry("CANCER_TYPE_2", 0.8, 1, 0),
                createTestPredictionEntry("CANCER_TYPE_3", 0.0, 3, 0)
        );

        CUPPA_PREDICTIONS = new CuppaPredictions(predictionEntries);
        DB_ACCESS.writeCuppa(TEST_SAMPLE_ID, CUPPA_PREDICTIONS, TOP_N_PROBS);

        CUPPA_RECORDS = DB_ACCESS.context()
                .select().from(Tables.CUPPA).fetch()
                .stream().map(r -> r.into(CuppaRecord.class))
                .toList();
    }

    @Test
    public void onlyTopNPredictionsAreLoaded()
    {
        assertEquals(6, CUPPA_PREDICTIONS.size());
        assertEquals(4, CUPPA_RECORDS.size());
    }

    @Test
    public void hasCorrectTopPredictions()
    {
        List<CuppaRecord> topPredictions = CUPPA_RECORDS.stream().filter(x -> x.getRank() == 1).toList();

        assertEquals("CANCER_TYPE_2", topPredictions.get(0).getCancertype());
        assertEquals(0.9, topPredictions.get(0).getProb(), 0.01);

        assertEquals("CANCER_TYPE_2", topPredictions.get(1).getCancertype());
        assertEquals(0.8, topPredictions.get(1).getProb(), 0.01);
    }

    private static CuppaPredictionEntry createTestPredictionEntry(String cancerType, double probability, int rank, int rankGroup)
    {
        return new CuppaPredictionEntry(
                TEST_SAMPLE_ID, DataType.PROB, ClassifierGroup.DNA, ClassifierName.DNA_COMBINED, null, Double.NaN,
                cancerType, probability, rank, rankGroup
        );
    }
}
