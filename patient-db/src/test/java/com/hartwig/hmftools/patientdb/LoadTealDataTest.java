package com.hartwig.hmftools.patientdb;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.teal.ImmutableTelomereLength;
import com.hartwig.hmftools.common.teal.TelomereLength;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.TelomerelengthRecord;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class LoadTealDataTest extends DatabaseTestBase
{
    @Test
    public void canWriteTelomereLengthData()
    {
        TelomereLength somaticTelomereLength = ImmutableTelomereLength.builder()
                .type("ref")
                .rawTelomereLength(3889.97)
                .finalTelomereLength(3889.97)
                .fullFragments(20433)
                .cRichPartialFragments(3396)
                .gRichPartialFragments(1729)
                .totalTelomericReads(42533)
                .purity(1)
                .ploidy(2)
                .duplicateProportion(0.105)
                .meanReadDepth(32.93)
                .gc50ReadDepth(33.6)
                .build();

        databaseAccess.writeTelomereLength(TEST_SAMPLE_ID, null, somaticTelomereLength);

        List<TelomerelengthRecord> telomereLengthRecords = fetchTable(Tables.TELOMERELENGTH, TelomerelengthRecord.class);
        assertEquals(1, telomereLengthRecords.size());
    }
}
