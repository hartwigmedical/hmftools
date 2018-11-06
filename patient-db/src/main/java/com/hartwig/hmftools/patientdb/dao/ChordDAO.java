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

        Double BRCA1 = chordAnalysis.BRCA1Value();
        Double none = chordAnalysis.noneValue();
        Double BRCA2 = chordAnalysis.BRCA2Value();
        Double hrd = chordAnalysis.hrdValue();
        Double predictedResponseValue = chordAnalysis.predictedResponseValue();

        context.insertInto(
                CHORD,
                CHORD.SAMPLEID,
                CHORD.BRCA1VALUE,
                CHORD.NONEVALUE,
                CHORD.BRCA2VALUE,
                CHORD.HRDVALUE,
                CHORD.PREDICTEDRESPONSEVALUE)
                .values(sample,
                        DatabaseUtil.decimal(BRCA1),
                        DatabaseUtil.decimal(none),
                        DatabaseUtil.decimal(BRCA2),
                        DatabaseUtil.decimal(hrd),
                        DatabaseUtil.decimal(predictedResponseValue))
                .execute();
    }

    void deleteChordForSample(@NotNull String sample) {
        context.delete(CHORD).where(CHORD.SAMPLEID.eq(sample)).execute();
    }
}
