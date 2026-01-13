package com.hartwig.hmftools.patientdb;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.chord.ImmutableChordData;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.ChordRecord;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class LoadChordDataTest extends DatabaseTestBase
{
    @Test
    public void canLoadChordData()
    {
        ChordData chordData = ImmutableChordData.builder()
                .BRCA1Value(0.25)
                .BRCA2Value(0.35)
                .hrdValue(0.60)
                .hrStatus(ChordStatus.HR_DEFICIENT)
                .hrdType("BRCA2_type")
                .remarksHrStatus("")
                .remarksHrdType("")
                .build();

        databaseAccess.writeChord("example", chordData);

        List<ChordRecord> chordRecords = fetchTable(Tables.CHORD, ChordRecord.class);

        assertEquals(1, chordRecords.size());

        ChordRecord actualRecord = chordRecords.get(0);
        ChordRecord expectedRecord = from(chordData, "example");

        assertEquals(expectedRecord, actualRecord);
    }

    private static ChordRecord from(final ChordData chordData, final String sampleId)
    {
        return new ChordRecord(
                sampleId,
                chordData.BRCA1Value(),
                chordData.BRCA2Value(),
                chordData.hrdValue(),
                chordData.hrStatus().toString(),
                chordData.hrdType(),
                chordData.remarksHrStatus(),
                chordData.remarksHrdType()
        );
    }
}
