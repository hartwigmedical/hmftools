package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PURITY;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PURITYRANGE;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;

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
                        .score(result.getValue(PURITY.SCORE))
                        .diploidProportion(result.getValue(PURITY.DIPLOIDPROPORTION))
                        .ploidy(result.getValue(PURITY.PLOIDY))
                        .build();
    }

    @Nullable
    PurityContext readPurityContext(@NotNull final String sample) {
        @Nullable
        Record result = context.select().from(PURITY).where(PURITY.SAMPLEID.eq(sample)).fetchOne();
        if (result == null) {
            return null;
        }

        final FittedPurity purity = ImmutableFittedPurity.builder()
                .purity(result.getValue(PURITY.PURITY_))
                .normFactor(result.getValue(PURITY.NORMFACTOR))
                .score(result.getValue(PURITY.SCORE))
                .diploidProportion(result.getValue(PURITY.DIPLOIDPROPORTION))
                .ploidy(result.getValue(PURITY.PLOIDY))
                .build();

        final FittedPurityScore score = ImmutableFittedPurityScore.builder()
                .minPurity(result.getValue(PURITY.MINPURITY))
                .maxPurity(result.getValue(PURITY.MAXPURITY))
                .minPloidy(result.getValue(PURITY.MINPLOIDY))
                .maxPloidy(result.getValue(PURITY.MAXPLOIDY))
                .minDiploidProportion(result.getValue(PURITY.MINDIPLOIDPROPORTION))
                .maxDiploidProportion(result.getValue(PURITY.MAXDIPLOIDPROPORTION))
                .build();

        return ImmutablePurityContext.builder()
                .bestFit(purity)
                .score(score)
                .gender(Gender.valueOf(result.getValue(PURITY.GENDER)))
                .polyClonalProportion(result.getValue(PURITY.POLYCLONALPROPORTION))
                .status(FittedPurityStatus.valueOf(result.getValue(PURITY.STATUS)))
                .build();
    }

    void write(@NotNull final String sample, @NotNull final PurityContext purity, @NotNull final PurpleQC checks) {

        final FittedPurity bestFit = purity.bestFit();
        final FittedPurityScore score = purity.score();

        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(PURITY).where(PURITY.SAMPLEID.eq(sample)).execute();

        context.insertInto(PURITY,
                PURITY.SAMPLEID,
                PURITY.PURITY_,
                PURITY.GENDER,
                PURITY.STATUS,
                PURITY.QCSTATUS,
                PURITY.NORMFACTOR,
                PURITY.SCORE,
                PURITY.PLOIDY,
                PURITY.DIPLOIDPROPORTION,
                PURITY.MINDIPLOIDPROPORTION,
                PURITY.MAXDIPLOIDPROPORTION,
                PURITY.MINPURITY,
                PURITY.MAXPURITY,
                PURITY.MINPLOIDY,
                PURITY.MAXPLOIDY,
                PURITY.POLYCLONALPROPORTION,
                PURITY.MODIFIED)
                .values(sample,
                        bestFit.purity(),
                        purity.gender().toString(),
                        purity.status().toString(),
                        checks.status().toString(),
                        bestFit.normFactor(),
                        bestFit.score(),
                        bestFit.ploidy(),
                        bestFit.diploidProportion(),
                        score.minDiploidProportion(),
                        score.maxDiploidProportion(),
                        score.minPurity(),
                        score.maxPurity(),
                        score.minPloidy(),
                        score.maxPloidy(),
                        purity.polyClonalProportion(),
                        timestamp)
                .execute();

    }

    void write(@NotNull final String sample, @NotNull List<FittedPurity> purities) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(PURITYRANGE).where(PURITYRANGE.SAMPLEID.eq(sample)).execute();

        InsertValuesStep7 inserter = context.insertInto(PURITYRANGE,
                PURITYRANGE.SAMPLEID,
                PURITYRANGE.PURITY,
                PURITYRANGE.NORMFACTOR,
                PURITYRANGE.SCORE,
                PURITYRANGE.PLOIDY,
                PURITYRANGE.DIPLOIDPROPORTION,
                PURITYRANGE.MODIFIED);

        purities.forEach(x -> addPurity(timestamp, inserter, sample, x));
        inserter.execute();
    }

    private void addPurity(Timestamp timestamp, InsertValuesStep7 inserter, String sample, FittedPurity purity) {
        inserter.values(sample,
                purity.purity(),
                purity.normFactor(),
                purity.score(),
                purity.ploidy(),
                purity.diploidProportion(),
                timestamp);
    }
}
