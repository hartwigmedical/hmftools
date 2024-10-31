package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CDR3LOCUSSUMMARY;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CDR3SEQUENCE;

import java.sql.Timestamp;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.common.cider.Cdr3LocusSummary;
import com.hartwig.hmftools.common.cider.Cdr3Sequence;
import com.hartwig.hmftools.common.cider.ImmutableCdr3LocusSummary;
import com.hartwig.hmftools.common.cider.ImmutableCdr3Sequence;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.Cdr3locussummaryRecord;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.Cdr3sequenceRecord;

import org.jooq.DSLContext;
import org.jooq.InsertValuesStep11;
import org.jooq.InsertValuesStep8;
import org.jooq.Record;
import org.jooq.Result;

public class CiderDAO
{
    private final DSLContext context;

    private static final int DB_BATCH_INSERT_SIZE = 10000;

    public CiderDAO(final DSLContext context) {
        this.context = context;
    }

    void deleteCiderDataForSample(final String sampleId)
    {
        context.delete(CDR3SEQUENCE).where(CDR3SEQUENCE.SAMPLEID.eq(sampleId)).execute();
        context.delete(CDR3LOCUSSUMMARY).where(CDR3LOCUSSUMMARY.SAMPLEID.eq(sampleId)).execute();
    }

    public void writeCdr3Sequence(final String sampleId, final Collection<Cdr3Sequence> cdr3Sequences)
    {
        context.delete(CDR3SEQUENCE).where(CDR3SEQUENCE.SAMPLEID.eq(sampleId)).execute();

        Timestamp timestamp = new Timestamp(new Date().getTime());

        try (InsertValuesStep11<Cdr3sequenceRecord, LocalDateTime, String, String, String, String, String, String, Integer, Integer, Byte, Byte> inserter
                = context.insertInto(CDR3SEQUENCE,
                CDR3SEQUENCE.MODIFIED,
                CDR3SEQUENCE.SAMPLEID,
                CDR3SEQUENCE.CDR3SEQ,
                CDR3SEQUENCE.CDR3AA,
                CDR3SEQUENCE.LOCUS,
                CDR3SEQUENCE.FILTER,
                CDR3SEQUENCE.BLASTNSTATUS,
                CDR3SEQUENCE.MINHIGHQUALBASEREADS,
                CDR3SEQUENCE.ASSIGNEDREADS,
                CDR3SEQUENCE.INFRAME,
                CDR3SEQUENCE.CONTAINSSTOP))
        {
            cdr3Sequences.forEach(cdr3Sequence ->
                inserter.values(timestamp.toLocalDateTime(),
                    sampleId,
                    DatabaseUtil.checkStringLength(cdr3Sequence.cdr3Seq(), CDR3SEQUENCE.CDR3SEQ),
                    DatabaseUtil.checkStringLength(cdr3Sequence.cdr3AA(), CDR3SEQUENCE.CDR3AA),
                    DatabaseUtil.checkStringLength(cdr3Sequence.locus(), CDR3SEQUENCE.LOCUS),
                    DatabaseUtil.checkStringLength(cdr3Sequence.filter(), CDR3SEQUENCE.FILTER),
                    DatabaseUtil.checkStringLength(cdr3Sequence.blastnStatus(), CDR3SEQUENCE.BLASTNSTATUS),
                    cdr3Sequence.minHighQualBaseReads(),
                    cdr3Sequence.assignedReads(),
                    DatabaseUtil.booleanToByte(cdr3Sequence.inFrame()),
                    DatabaseUtil.booleanToByte(cdr3Sequence.containsStop())));

            inserter.execute();
        }
    }

    public List<Cdr3Sequence> readCdr3Sequence(final String sampleId)
    {
        List<Cdr3Sequence> result = new ArrayList<>();
        Result<Record> queryResult = context.select().from(CDR3SEQUENCE).where(CDR3SEQUENCE.SAMPLEID.eq(sampleId)).fetch();

        for(Record record : queryResult)
        {
            result.add(ImmutableCdr3Sequence.builder()
                    .cdr3Seq(record.get(CDR3SEQUENCE.CDR3SEQ))
                    .cdr3AA(record.get(CDR3SEQUENCE.CDR3AA))
                    .locus(record.get(CDR3SEQUENCE.LOCUS))
                    .filter(record.get(CDR3SEQUENCE.FILTER))
                    .blastnStatus(record.get(CDR3SEQUENCE.BLASTNSTATUS))
                    .minHighQualBaseReads(record.get(CDR3SEQUENCE.MINHIGHQUALBASEREADS))
                    .assignedReads(record.get(CDR3SEQUENCE.ASSIGNEDREADS))
                    .inFrame(DatabaseUtil.byteToBoolean(record.get(CDR3SEQUENCE.INFRAME)))
                    .containsStop(DatabaseUtil.byteToBoolean(record.get(CDR3SEQUENCE.CONTAINSSTOP)))
                    .build());
        }
        return result;
    }

    public void writeLocusSummaries(final String sampleId, final Collection<Cdr3LocusSummary> locusSummaries)
    {
        context.delete(CDR3LOCUSSUMMARY).where(CDR3LOCUSSUMMARY.SAMPLEID.eq(sampleId)).execute();

        Timestamp timestamp = new Timestamp(new Date().getTime());

        try (InsertValuesStep8<Cdr3locussummaryRecord, LocalDateTime, String, String, Integer, Integer, Byte, Integer, Integer> inserter
                = context.insertInto(CDR3LOCUSSUMMARY,
                        CDR3LOCUSSUMMARY.MODIFIED,
                        CDR3LOCUSSUMMARY.SAMPLEID,
                        CDR3LOCUSSUMMARY.LOCUS,
                        CDR3LOCUSSUMMARY.READSUSED,
                        CDR3LOCUSSUMMARY.READSTOTAL,
                        CDR3LOCUSSUMMARY.DOWNSAMPLED,
                        CDR3LOCUSSUMMARY.SEQUENCES,
                        CDR3LOCUSSUMMARY.PASSSEQUENCES))
        {
                locusSummaries.forEach( o -> inserter.values(timestamp.toLocalDateTime(),
                            sampleId,
                            DatabaseUtil.checkStringLength(o.locus(), CDR3LOCUSSUMMARY.LOCUS),
                            o.readsUsed(),
                            o.readsTotal(),
                            DatabaseUtil.booleanToByte(o.downSampled()),
                            o.sequences(),
                            o.passSequences()));

                inserter.execute();
        }
    }

    public List<Cdr3LocusSummary> readCdr3LocusSummaries(final String sampleId)
    {
        List<Cdr3LocusSummary> result = new ArrayList<>();
        Result<Record> queryResult = context.select().from(CDR3LOCUSSUMMARY).where(CDR3LOCUSSUMMARY.SAMPLEID.eq(sampleId)).fetch();

        for(Record record : queryResult)
        {
            result.add(ImmutableCdr3LocusSummary.builder()
                            .locus(record.get(CDR3LOCUSSUMMARY.LOCUS))
                            .readsUsed(record.get(CDR3LOCUSSUMMARY.READSUSED))
                            .readsTotal(record.get(CDR3LOCUSSUMMARY.READSTOTAL))
                            .downSampled(DatabaseUtil.byteToBoolean(record.get(CDR3LOCUSSUMMARY.DOWNSAMPLED)))
                            .sequences(record.get(CDR3LOCUSSUMMARY.SEQUENCES))
                            .passSequences(record.get(CDR3LOCUSSUMMARY.PASSSEQUENCES))
                    .build());
        }
        return result;
    }
}
