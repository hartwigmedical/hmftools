package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CHORD;

import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordStatus;
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
                        DatabaseUtil.decimal(chordAnalysis.BRCA1Value()),
                        DatabaseUtil.decimal(chordAnalysis.BRCA2Value()),
                        DatabaseUtil.decimal(chordAnalysis.hrdValue()),
                        chordAnalysis.hrStatus().toString(),
                        chordAnalysis.hrdType(),
                        chordAnalysis.remarksHrStatus(),
                        chordAnalysis.remarksHrdType())
                .execute();
    }

    @Nullable
    ChordAnalysis readChordAnalysis(@NotNull String sample) {
        Record result = context.select().from(CHORD).where(CHORD.SAMPLEID.eq(sample)).fetchOne();

        if (result == null) {
            return null;
        }

        return ImmutableChordAnalysis.builder()
                .BRCA1Value(result.getValue(CHORD.BRCA1))
                .BRCA2Value(result.getValue(CHORD.BRCA2))
                .hrdValue(result.getValue(CHORD.HRD))
                .hrStatus(ChordStatus.valueOf(result.getValue(CHORD.HRSTATUS)))
                .hrdType(result.getValue(CHORD.HRDTYPE))
                .remarksHrStatus(result.getValue(CHORD.REMARKSHRSTATUS))
                .remarksHrdType(result.getValue(CHORD.REMARKSHRDTYPE))
                .build();
    }

    void deleteChordForSample(@NotNull String sample) {
        context.delete(CHORD).where(CHORD.SAMPLEID.eq(sample)).execute();
    }
}
