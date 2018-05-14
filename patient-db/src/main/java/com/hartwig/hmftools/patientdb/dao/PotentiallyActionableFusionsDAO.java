package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.BASELINE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SAMPLE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTBREAKEND;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTFUSION;

import java.sql.ResultSet;
import java.util.stream.Stream;

import com.hartwig.hmftools.patientdb.data.PotentialActionableFusion;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Structuralvariantbreakend;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

public class PotentiallyActionableFusionsDAO {
    public static final Structuralvariantbreakend FIVE_BREAKEND = STRUCTURALVARIANTBREAKEND.as("fiveBreakend");
    public static final Structuralvariantbreakend THREE_BREAKEND = STRUCTURALVARIANTBREAKEND.as("threeBreakend");

    @NotNull
    private final DSLContext context;

    PotentiallyActionableFusionsDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @NotNull
    Stream<PotentialActionableFusion> potentiallyActionableFusions() {
        return context.select(STRUCTURALVARIANT.SAMPLEID, BASELINE.PRIMARYTUMORLOCATION, FIVE_BREAKEND.GENE, THREE_BREAKEND.GENE)
                .from(STRUCTURALVARIANTFUSION.join(FIVE_BREAKEND)
                        .on(FIVE_BREAKEND.ID.eq(STRUCTURALVARIANTFUSION.FIVEPRIMEBREAKENDID))
                        .join(THREE_BREAKEND)
                        .on(THREE_BREAKEND.ID.eq(STRUCTURALVARIANTFUSION.THREEPRIMEBREAKENDID))
                        .join(STRUCTURALVARIANT)
                        .on(STRUCTURALVARIANT.ID.eq(FIVE_BREAKEND.STRUCTURALVARIANTID))
                        .join(SAMPLE)
                        .on(STRUCTURALVARIANT.SAMPLEID.eq(SAMPLE.SAMPLEID))
                        .join(BASELINE)
                        .on(SAMPLE.PATIENTID.eq(BASELINE.PATIENTID)))
                .where(STRUCTURALVARIANTFUSION.ISREPORTED.eq((byte) 1))
                .fetchSize(Integer.MIN_VALUE)
                .resultSetType(ResultSet.TYPE_FORWARD_ONLY)
                .resultSetConcurrency(ResultSet.CONCUR_READ_ONLY)
                .stream()
                .map(PotentialActionableFusion::of);
    }
}
