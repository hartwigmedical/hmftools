package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CHORD;

import com.hartwig.hmftools.common.chord.ChordAnalysis;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

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

    void deleteChordForSample(@NotNull String sample) {
        context.delete(CHORD).where(CHORD.SAMPLEID.eq(sample)).execute();
    }
}
