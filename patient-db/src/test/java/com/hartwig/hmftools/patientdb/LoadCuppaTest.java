package com.hartwig.hmftools.patientdb;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictions;
import java.io.IOException;

public class LoadCuppaTest {
    /*
    mysql --user="writer"  --password

    USE patientdb;
    source /Users/lnguyen/Hartwig/hartwigmedical/hmftools/patient-db/src/main/resources/generate_database.sql;
    */

    private static final String CUPPA_VIS_DATA_TSV = "/Users/lnguyen/Hartwig/hartwigmedical/hmftools/hmf-common/src/test/resources/cuppa/cuppa_vis_data.tsv";

    public static void main(String[] args) throws IOException
    {
        CuppaPredictions cuppaPredictions = CuppaPredictions.fromTsv(CUPPA_VIS_DATA_TSV);
        cuppaPredictions.getTopPredictions(3).sortByRank().printPredictions();
    }
}
