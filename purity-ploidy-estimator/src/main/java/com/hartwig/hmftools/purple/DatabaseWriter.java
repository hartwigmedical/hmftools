package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PURITY;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Copynumber.COPYNUMBER;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityScore;
import com.hartwig.hmftools.common.purple.region.ConsolidatedRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep9;
import org.jooq.SQLDialect;
import org.jooq.impl.DSL;

class DatabaseWriter {

    private static final Logger LOGGER = LogManager.getLogger(DatabaseWriter.class);

    @NotNull
    private final DSLContext context;
    @NotNull
    private final Timestamp timestamp = new Timestamp(new Date().getTime());

    DatabaseWriter(@NotNull final String userName, @NotNull final String password, @NotNull final String url)
            throws SQLException {
        final Connection conn = DriverManager.getConnection(url, userName, password);
        this.context = DSL.using(conn, SQLDialect.MYSQL);
    }

    void writePurity(@NotNull final String sample, @NotNull FittedPurity purity, @NotNull FittedPurityScore score) {
        context.delete(PURITY).where(PURITY.SAMPLEID.eq(sample)).execute();

        context.insertInto(PURITY, PURITY.SAMPLEID, PURITY.PURITY_, PURITY.NORMFACTOR, PURITY.SCORE, PURITY.PLOIDY,
                PURITY.DIPLOIDPROPORTION, PURITY.POLYCLONALPROPORTION, PURITY.MINPURITY, PURITY.MAXPURITY,
                PURITY.MINPLOIDY, PURITY.MAXPLOIDY, PURITY.MODIFIED)
                .values(sample, purity.purity(), purity.normFactor(), purity.score(), purity.ploidy(),
                        purity.diplodProportion(), score.polyclonalProportion(), score.minPurity(), score.maxPurity(),
                        score.minPloidy(), score.maxPloidy(), timestamp)
                .execute();

    }

    void writeCopynumbers(@NotNull final String sample, @NotNull List<ConsolidatedRegion> regions) {
        context.delete(COPYNUMBER).where(COPYNUMBER.SAMPLEID.eq(sample)).execute();

        InsertValuesStep9 inserter = context.insertInto(COPYNUMBER, COPYNUMBER.SAMPLEID, COPYNUMBER.CHROMOSOME,
                COPYNUMBER.START, COPYNUMBER.END, COPYNUMBER.BAFCOUNT, COPYNUMBER.OBSERVEDBAF, COPYNUMBER.ACTUALBAF,
                COPYNUMBER.COPYNUMBER_, COPYNUMBER.MODIFIED);

        regions.forEach(x -> addCopynumberRecord(inserter, sample, x));
        inserter.execute();
    }

    private void addCopynumberRecord(InsertValuesStep9 inserter, String sample, ConsolidatedRegion region) {
        inserter.values(sample, region.chromosome(), region.start(), region.end(), region.bafCount(),
                region.averageObservedBAF(), region.averagePurityAdjustedBAF(), region.averageTumorCopyNumber(),
                timestamp);
    }
}
