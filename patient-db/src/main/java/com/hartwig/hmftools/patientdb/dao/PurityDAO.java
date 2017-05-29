package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PURITY;

import java.sql.Timestamp;
import java.util.Date;

import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityScore;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.Record;

class PurityDAO {

    @NotNull
    private final DSLContext context;

    PurityDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @Nullable
    public FittedPurity readFittedPurity(@NotNull final String sample) {
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

    public void write(@NotNull final String sample, @NotNull FittedPurity purity, @NotNull FittedPurityScore score) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(PURITY).where(PURITY.SAMPLEID.eq(sample)).execute();

        context.insertInto(PURITY, PURITY.SAMPLEID, PURITY.PURITY_, PURITY.NORMFACTOR, PURITY.SCORE, PURITY.PLOIDY,
                PURITY.DIPLOIDPROPORTION, PURITY.POLYCLONALPROPORTION, PURITY.MINPURITY, PURITY.MAXPURITY,
                PURITY.MINPLOIDY, PURITY.MAXPLOIDY, PURITY.MODIFIED)
                .values(sample, purity.purity(), purity.normFactor(), purity.score(), purity.ploidy(),
                        purity.diploidProportion(), score.polyclonalProportion(), score.minPurity(), score.maxPurity(),
                        score.minPloidy(), score.maxPloidy(), timestamp)
                .execute();

    }
}
