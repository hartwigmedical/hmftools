package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.BASELINE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GENECOPYNUMBER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PURITY;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SAMPLE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTBREAKEND;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTFUSION;

import java.sql.ResultSet;
import java.util.Collection;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.copynumber.CopyNumberAlteration;
import com.hartwig.hmftools.patientdb.data.PotentialActionableCNV;
import com.hartwig.hmftools.patientdb.data.PotentialActionableFusion;
import com.hartwig.hmftools.patientdb.data.PotentialActionableVariant;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Structuralvariantbreakend;

import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.ResultQuery;
import org.jooq.SelectConditionStep;

public class PotentiallyActionableItemsDAO {
    public static final Structuralvariantbreakend FIVE_BREAKEND = STRUCTURALVARIANTBREAKEND.as("fiveBreakend");
    public static final Structuralvariantbreakend THREE_BREAKEND = STRUCTURALVARIANTBREAKEND.as("threeBreakend");

    @NotNull
    private final DSLContext context;

    PotentiallyActionableItemsDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @NotNull
    Stream<PotentialActionableVariant> potentiallyActionableVariants(@NotNull final Collection<String> samples) {
        final ResultQuery<?> query = potentiallyActionableVariantsQuery(samples);
        return streamResults(query).map(PotentialActionableVariant::of);
    }

    @NotNull
    Stream<PotentialActionableCNV> potentiallyActionableCNVs(@NotNull final Collection<String> samples) {
        final ResultQuery<?> query = potentiallyActionableCNVsQuery(samples);
        return streamResults(query).map(PotentialActionableCNV::of).filter(cnv -> cnv.alteration() != CopyNumberAlteration.NEUTRAL);
    }

    @NotNull
    Stream<PotentialActionableFusion> potentiallyActionableFusions(@NotNull final Collection<String> samples) {
        final ResultQuery<?> query = potentiallyActionableFusionsQuery(samples);
        return streamResults(query).map(PotentialActionableFusion::of);
    }

    @NotNull
    private ResultQuery<?> potentiallyActionableVariantsQuery(@NotNull final Collection<String> samples) {
        SelectConditionStep<?> query = context.select(SOMATICVARIANT.SAMPLEID,
                SOMATICVARIANT.GENE,
                SOMATICVARIANT.CHROMOSOME,
                SOMATICVARIANT.POSITION,
                SOMATICVARIANT.REF,
                SOMATICVARIANT.ALT,
                SOMATICVARIANT.WORSTEFFECTTRANSCRIPT,
                SOMATICVARIANT.WORSTEFFECT,
                SOMATICVARIANT.TYPE)
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.FILTER.eq("PASS")
                        .and(SOMATICVARIANT.WORSTCODINGEFFECT.ne("NONE"))
                        .and(SOMATICVARIANT.WORSTCODINGEFFECT.ne("SYNONYMOUS")));
        if (samples.size() > 0) {
            query.and(SOMATICVARIANT.SAMPLEID.in(samples));
        }
        return query;
    }

    @NotNull
    private ResultQuery<?> potentiallyActionableCNVsQuery(@NotNull final Collection<String> samples) {
        SelectConditionStep<?> query = context.select(GENECOPYNUMBER.SAMPLEID, GENECOPYNUMBER.GENE, GENECOPYNUMBER.MINCOPYNUMBER)
                .from(GENECOPYNUMBER.join(PURITY).on(GENECOPYNUMBER.SAMPLEID.eq(PURITY.SAMPLEID)))
                .where(PURITY.QCSTATUS.eq("PASS"))
                .and(PURITY.STATUS.ne("NO_TUMOR"))
                .and(PURITY.PURITY_.ge(0.2))
                .and(GENECOPYNUMBER.MINCOPYNUMBER.le(0.5).or(GENECOPYNUMBER.MINCOPYNUMBER.div(PURITY.PLOIDY).ge(3.0)));
        if (samples.size() > 0) {
            query.and(GENECOPYNUMBER.SAMPLEID.in(samples));
        }
        return query;
    }

    @NotNull
    private ResultQuery<?> potentiallyActionableFusionsQuery(@NotNull final Collection<String> samples) {
        SelectConditionStep<?> query = context.select(STRUCTURALVARIANT.SAMPLEID, FIVE_BREAKEND.GENE, THREE_BREAKEND.GENE)
                .from(STRUCTURALVARIANTFUSION.join(FIVE_BREAKEND)
                        .on(FIVE_BREAKEND.ID.eq(STRUCTURALVARIANTFUSION.FIVEPRIMEBREAKENDID))
                        .join(THREE_BREAKEND)
                        .on(THREE_BREAKEND.ID.eq(STRUCTURALVARIANTFUSION.THREEPRIMEBREAKENDID))
                        .join(STRUCTURALVARIANT)
                        .on(STRUCTURALVARIANT.ID.eq(FIVE_BREAKEND.STRUCTURALVARIANTID)))
                .where(STRUCTURALVARIANTFUSION.ISREPORTED.eq((byte) 1));
        if (samples.size() > 0) {
            query.and(STRUCTURALVARIANT.SAMPLEID.in(samples));
        }
        return query;
    }

    @NotNull
    Stream<Pair<String, String>> allSampleAndTumorLocations(@NotNull final String sampleId) {
        return context.select(SAMPLE.SAMPLEID, BASELINE.PRIMARYTUMORLOCATION)
                .from(SAMPLE.leftJoin(BASELINE).on(SAMPLE.PATIENTID.eq(BASELINE.PATIENTID)))
                .where(SAMPLE.SAMPLEID.eq(sampleId))
                .fetch()
                .stream()
                .map(result -> Pair.create(result.get(SAMPLE.SAMPLEID), result.get(BASELINE.PRIMARYTUMORLOCATION)))
                .distinct();

    }

    private static <T extends Record> Stream<T> streamResults(@NotNull final ResultQuery<T> query) {
        return query.fetchSize(Integer.MIN_VALUE)
                .resultSetType(ResultSet.TYPE_FORWARD_ONLY)
                .resultSetConcurrency(ResultSet.CONCUR_READ_ONLY)
                .stream();
    }
}
