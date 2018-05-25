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
import java.util.stream.Stream;

import com.hartwig.hmftools.common.copynumber.CopyNumberAlteration;
import com.hartwig.hmftools.patientdb.data.PotentialActionableCNV;
import com.hartwig.hmftools.patientdb.data.PotentialActionableFusion;
import com.hartwig.hmftools.patientdb.data.PotentialActionableVariant;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Structuralvariantbreakend;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

public class PotentiallyActionableItemsDAO {
    public static final Structuralvariantbreakend FIVE_BREAKEND = STRUCTURALVARIANTBREAKEND.as("fiveBreakend");
    public static final Structuralvariantbreakend THREE_BREAKEND = STRUCTURALVARIANTBREAKEND.as("threeBreakend");

    @NotNull
    private final DSLContext context;

    PotentiallyActionableItemsDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @NotNull
    Stream<PotentialActionableVariant> potentiallyActionableVariants() {
        return context.select(SOMATICVARIANT.SAMPLEID,
                BASELINE.PRIMARYTUMORLOCATION,
                SOMATICVARIANT.GENE,
                SOMATICVARIANT.CHROMOSOME,
                SOMATICVARIANT.POSITION,
                SOMATICVARIANT.REF,
                SOMATICVARIANT.ALT)
                .from(SOMATICVARIANT.join(SAMPLE)
                        .on(SOMATICVARIANT.SAMPLEID.eq(SAMPLE.SAMPLEID))
                        .leftJoin(BASELINE)
                        .on(SAMPLE.PATIENTID.eq(BASELINE.PATIENTID)))
                .where(SOMATICVARIANT.FILTER.eq("PASS"))
                .fetchSize(Integer.MIN_VALUE)
                .resultSetType(ResultSet.TYPE_FORWARD_ONLY)
                .resultSetConcurrency(ResultSet.CONCUR_READ_ONLY)
                .stream()
                .map(PotentialActionableVariant::of);
    }

    @NotNull
    Stream<PotentialActionableCNV> potentiallyActionableCNVs() {
        return context.select(GENECOPYNUMBER.SAMPLEID, BASELINE.PRIMARYTUMORLOCATION, GENECOPYNUMBER.GENE, GENECOPYNUMBER.MINCOPYNUMBER)
                .from(GENECOPYNUMBER.join(PURITY)
                        .on(GENECOPYNUMBER.SAMPLEID.eq(PURITY.SAMPLEID))
                        .join(SAMPLE)
                        .on(GENECOPYNUMBER.SAMPLEID.eq(SAMPLE.SAMPLEID))
                        .leftJoin(BASELINE)
                        .on(SAMPLE.PATIENTID.eq(BASELINE.PATIENTID)))
                .where(PURITY.QCSTATUS.eq("PASS")
                        .and(PURITY.STATUS.ne("NO_TUMOR"))
                        .and(PURITY.PURITY_.ge(0.2))
                        .and(GENECOPYNUMBER.MINCOPYNUMBER.le(0.5).or(GENECOPYNUMBER.MINCOPYNUMBER.div(PURITY.PLOIDY).ge(3.0))))
                .fetchSize(Integer.MIN_VALUE)
                .resultSetType(ResultSet.TYPE_FORWARD_ONLY)
                .resultSetConcurrency(ResultSet.CONCUR_READ_ONLY)
                .stream()
                .map(PotentialActionableCNV::of)
                .filter(cnv -> cnv.alteration() != CopyNumberAlteration.NEUTRAL);
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
                        .leftJoin(BASELINE)
                        .on(SAMPLE.PATIENTID.eq(BASELINE.PATIENTID)))
                .where(STRUCTURALVARIANTFUSION.ISREPORTED.eq((byte) 1))
                .fetchSize(Integer.MIN_VALUE)
                .resultSetType(ResultSet.TYPE_FORWARD_ONLY)
                .resultSetConcurrency(ResultSet.CONCUR_READ_ONLY)
                .stream()
                .map(PotentialActionableFusion::of);
    }
}
