package com.hartwig.hmftools.patientdb;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.cider.Cdr3LocusSummary;
import com.hartwig.hmftools.common.cider.ImmutableCdr3LocusSummary;
import com.hartwig.hmftools.common.cider.ImmutableCdr3Sequence;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.Cdr3locussummaryRecord;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.Cdr3sequenceRecord;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class LoadCiderDataTest extends DatabaseTestBase
{
    private static final String TEST_SAMPLE_ID = "example";

    @Test
    public void canWriteCdr3Sequence()
    {
        ImmutableCdr3Sequence sequence = ImmutableCdr3Sequence.builder()
                .cdr3Seq("CCT")
                .cdr3AA("P")
                .locus("IGL")
                .filter("NO_VDJ_ALIGNMENT;MIN_LENGTH")
                .blastnStatus("NO_VDJ_ALIGNMENT")
                .minHighQualBaseReads(1)
                .assignedReads(2)
                .inFrame(true)
                .containsStop(false)
                .build();

        databaseAccess.writeCdr3Sequences(TEST_SAMPLE_ID, List.of(sequence));
        List<Cdr3sequenceRecord> records = fetchTable(Tables.CDR3SEQUENCE, Cdr3sequenceRecord.class);
        assertEquals(1, records.size());
    }

    @Test
    public void canWriteCdr3LocusSummary()
    {
        Cdr3LocusSummary summary = ImmutableCdr3LocusSummary.builder()
                .locus("IGH")
                .readsUsed(18221)
                .readsTotal(18685)
                .downSampled(false)
                .sequences(29)
                .passSequences(0)
                .build();

        databaseAccess.writeCdr3LocusSummaries(TEST_SAMPLE_ID, List.of(summary));
        List<Cdr3locussummaryRecord> records = fetchTable(Tables.CDR3LOCUSSUMMARY, Cdr3locussummaryRecord.class);
        assertEquals(1, records.size());
    }
}
