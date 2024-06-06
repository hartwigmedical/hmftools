package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TELOMERELENGTH;

import java.sql.Timestamp;
import java.time.LocalDateTime;
import java.util.Date;

import com.hartwig.hmftools.common.teal.TelomereLength;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.TelomerelengthRecord;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.DeleteUsingStep;
import org.jooq.InsertValuesStep11;

public class TealDAO
{
    private final DSLContext context;

    public TealDAO(final DSLContext context) {
        this.context = context;
    }

    void deleteTealDataForSample(final String sampleId)
    {
        try (DeleteUsingStep<TelomerelengthRecord> deleter = context.delete(TELOMERELENGTH))
        {
            deleter.where(TELOMERELENGTH.SAMPLEID.eq(sampleId)).execute();
        }
    }

    public void writeTelomereLength(final String sampleId, @Nullable TelomereLength germelineTelomereLength,
            @Nullable TelomereLength somaticTelomereLength)
    {
        deleteTealDataForSample(sampleId);

        Timestamp timestamp = new Timestamp(new Date().getTime());

        try (InsertValuesStep11<TelomerelengthRecord, LocalDateTime, String, Double, Double, Integer, Integer, Integer, Integer, Integer, Integer, Double> inserter
                = context.insertInto(TELOMERELENGTH,
                TELOMERELENGTH.MODIFIED,
                TELOMERELENGTH.SAMPLEID,
                TELOMERELENGTH.GERMLINETELOMERELENGTH,
                TELOMERELENGTH.SOMATICTELOMERELENGTH,
                TELOMERELENGTH.GERMLINEFULLFRAGMENTS,
                TELOMERELENGTH.GERMLINECRICHPARTIALFRAGMENTS,
                TELOMERELENGTH.GERMLINEGRICHPARTIALFRAGMENTS,
                TELOMERELENGTH.SOMATICFULLFRAGMENTS,
                TELOMERELENGTH.SOMATICCRICHPARTIALFRAGMENTS,
                TELOMERELENGTH.SOMATICGRICHPARTIALFRAGMENTS,
                TELOMERELENGTH.SAMPLEMIXLENGTH))
        {
            inserter.values(timestamp.toLocalDateTime(),
                    sampleId,
                    germelineTelomereLength != null ? germelineTelomereLength.finalTelomereLength() : null,
                    somaticTelomereLength != null ? somaticTelomereLength.finalTelomereLength() : null,
                    germelineTelomereLength != null ? germelineTelomereLength.fullFragments() : null,
                    germelineTelomereLength != null ? germelineTelomereLength.cRichPartialFragments() : null,
                    germelineTelomereLength != null ? germelineTelomereLength.gRichPartialFragments() : null,
                    somaticTelomereLength != null ? somaticTelomereLength.fullFragments() : null,
                    somaticTelomereLength != null ? somaticTelomereLength.cRichPartialFragments() : null,
                    somaticTelomereLength != null ? somaticTelomereLength.gRichPartialFragments() : null,
                    somaticTelomereLength != null ? somaticTelomereLength.rawTelomereLength() : null);

            inserter.execute();
        }
    }
}
