package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PURITY;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PURITYRANGE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PURITYSCORE;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep7;
import org.jooq.Record;

class PurityDAO {

    @NotNull
    private final DSLContext context;

    PurityDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @Nullable
    FittedPurity readFittedPurity(@NotNull final String sample) {
        @Nullable
        Record result = context.select().from(PURITY).where(PURITY.SAMPLEID.eq(sample)).fetchOne();
        return result == null
                ? null
                : ImmutableFittedPurity.builder()
                        .purity(result.getValue(PURITY.PURITY_))
                        .normFactor(result.getValue(PURITY.NORMFACTOR))
                        .modelBAFDeviation(0)
                        .score(result.getValue(PURITY.SCORE))
                        .diploidProportion(result.getValue(PURITY.DIPLOIDPROPORTION))
                        .build();
    }

    void write(@NotNull final String sample, @NotNull FittedPurityScore score) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(PURITYSCORE).where(PURITYSCORE.SAMPLEID.eq(sample)).execute();

        context.insertInto(PURITYSCORE, PURITYSCORE.SAMPLEID, PURITYSCORE.POLYCLONALPROPORTION, PURITYSCORE.MINPURITY,
                PURITYSCORE.MAXPURITY, PURITYSCORE.MINPLOIDY, PURITYSCORE.MAXPLOIDY, PURITYSCORE.MODIFIED)
                .values(sample, score.polyclonalProportion(), score.minPurity(), score.maxPurity(), score.minPloidy(),
                        score.maxPloidy(), timestamp)
                .execute();
    }

    void write(@NotNull final String sample, @NotNull List<FittedPurity> purities) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(PURITYRANGE).where(PURITYRANGE.SAMPLEID.eq(sample)).execute();

        InsertValuesStep7 inserter = context.insertInto(PURITYRANGE, PURITYRANGE.SAMPLEID, PURITYRANGE.PURITY, PURITYRANGE.NORMFACTOR,
                PURITYRANGE.SCORE, PURITYRANGE.PLOIDY, PURITYRANGE.DIPLOIDPROPORTION, PURITYRANGE.MODIFIED);

        purities.forEach(x -> addPurity(timestamp, inserter, sample, x));
        inserter.execute();
    }

    private void addPurity(Timestamp timestamp, InsertValuesStep7 inserter, String sample, FittedPurity purity) {
        inserter.values(sample, purity.purity(), purity.normFactor(), purity.score(), purity.ploidy(),
                purity.diploidProportion(), timestamp);
    }
}
