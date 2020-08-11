package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CHORD;

import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ImmutableChordAnalysis;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.Record;

class ChordDAO {

    @NotNull
    private final DSLContext context;

    ChordDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeChord(@NotNull String sample, @NotNull ChordAnalysis chordAnalysis) {
        deleteChordForSample(sample);

        context.insertInto(CHORD, CHORD.SAMPLEID, CHORD.NOTHING, CHORD.BRCA1, CHORD.BRCA2, CHORD.HRD, CHORD.PREDICTEDRESPONSE)
                .values(sample,
                        DatabaseUtil.decimal(chordAnalysis.noneValue()),
                        DatabaseUtil.decimal(chordAnalysis.BRCA1Value()),
                        DatabaseUtil.decimal(chordAnalysis.BRCA2Value()),
                        DatabaseUtil.decimal(chordAnalysis.hrdValue()),
                        chordAnalysis.predictedResponseValue() ? (byte) 1 : (byte) 0)
                .execute();
    }

    @Nullable
    ChordAnalysis readChordAnalysis(@NotNull String sample) {
        Record result = context.select().from(CHORD).where(CHORD.SAMPLEID.eq(sample)).fetchOne();

        if (result == null) {
            return null;
        }

        return ImmutableChordAnalysis.builder()
                .noneValue(result.getValue(CHORD.NOTHING))
                .BRCA1Value(result.getValue(CHORD.BRCA1))
                .BRCA2Value(result.getValue(CHORD.BRCA2))
                .hrdValue(result.getValue(CHORD.HRD))
                .predictedResponseValue(result.getValue(CHORD.PREDICTEDRESPONSE)==1)
                .build();
    }


    void deleteChordForSample(@NotNull String sample) {
        context.delete(CHORD).where(CHORD.SAMPLEID.eq(sample)).execute();
    }
}
