package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CHORD;

import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.chord.ImmutableChordData;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.Record;

class ChordDAO {

    private final DSLContext context;

    ChordDAO(final DSLContext context) {
        this.context = context;
    }

    public ChordData readChord(final String sampleId)
    {
        Record result = context.select().from(CHORD).where(CHORD.SAMPLEID.eq(sampleId)).fetchOne();
        if (result == null) {
            return null;
        }

        return ImmutableChordData.builder()
                .hrdValue(result.getValue(CHORD.HRD))
                .BRCA1Value(result.getValue(CHORD.BRCA1))
                .BRCA2Value(result.getValue(CHORD.BRCA2))
                .hrStatus(ChordStatus.valueOf(result.getValue(CHORD.HRSTATUS)))
                .hrdType(result.getValue(CHORD.HRDTYPE))
                .remarksHrdType(result.getValue(CHORD.REMARKSHRDTYPE))
                .remarksHrStatus(result.getValue(CHORD.REMARKSHRSTATUS))
                .build();
    }

    void writeChord(final String sample, final ChordData chordData)
    {
        deleteChordForSample(sample);

        context.insertInto(CHORD,
                CHORD.SAMPLEID,
                CHORD.BRCA1,
                CHORD.BRCA2,
                CHORD.HRD,
                CHORD.HRSTATUS,
                CHORD.HRDTYPE,
                CHORD.REMARKSHRSTATUS,
                CHORD.REMARKSHRDTYPE)
                .values(sample,
                        DatabaseUtil.decimal(chordData.BRCA1Value()),
                        DatabaseUtil.decimal(chordData.BRCA2Value()),
                        DatabaseUtil.decimal(chordData.hrdValue()),
                        chordData.hrStatus().toString(),
                        chordData.hrdType(),
                        chordData.remarksHrStatus(),
                        chordData.remarksHrdType())
                .execute();
    }

    void deleteChordForSample(final String sample) {
        context.delete(CHORD).where(CHORD.SAMPLEID.eq(sample)).execute();
    }
}
