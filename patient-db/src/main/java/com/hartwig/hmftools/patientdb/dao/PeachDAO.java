package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PEACHCALLS;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PEACHGENOTYPE;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.peach.PeachCalls;
import com.hartwig.hmftools.common.peach.PeachGenotype;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep12;
import org.jooq.InsertValuesStep9;

class PeachDAO {

    @NotNull
    private final DSLContext context;

    PeachDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writePeach(@NotNull String sample, @NotNull List<PeachGenotype> peachGenotypes, @NotNull List<PeachCalls> peachCalls) {
        deletePeachForSample(sample);

        Timestamp timestamp = new Timestamp(new Date().getTime());

        for (List<PeachGenotype> genotypes : Iterables.partition(peachGenotypes, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep9 inserter = context.insertInto(PEACHGENOTYPE,
                    PEACHGENOTYPE.MODIFIED,
                    PEACHGENOTYPE.SAMPLEID,
                    PEACHGENOTYPE.GENE,
                    PEACHGENOTYPE.HAPLOTYPE,
                    PEACHGENOTYPE.FUNCTION,
                    PEACHGENOTYPE.LINKEDDRUGS,
                    PEACHGENOTYPE.URLPRESCRIPTIONINFO,
                    PEACHGENOTYPE.PANELVERSION,
                    PEACHGENOTYPE.REPOVERSION);
            genotypes.forEach(x -> addGenotype(timestamp, inserter, sample, x));
            inserter.execute();
        }

        for (List<PeachCalls> calls : Iterables.partition(peachCalls, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep12 inserter = context.insertInto(PEACHCALLS,
                    PEACHCALLS.MODIFIED,
                    PEACHCALLS.SAMPLEID,
                    PEACHCALLS.GENE,
                    PEACHCALLS.POSITIONGRCH37,
                    PEACHCALLS.REFGRCH37,
                    PEACHCALLS.ALTGRCH37,
                    PEACHCALLS.POSITIONGRCH38,
                    PEACHCALLS.REFGRCH38,
                    PEACHCALLS.ALTGRCH38,
                    PEACHCALLS.RSID,
                    PEACHCALLS.VARIANTANNOTATION,
                    PEACHCALLS.FILTER);
            calls.forEach(x -> addCalls(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addGenotype(@NotNull Timestamp timestamp, @NotNull InsertValuesStep9 inserter, @NotNull String sample,
            @NotNull PeachGenotype genotype) {
        inserter.values(timestamp,
                sample,
                genotype.gene(),
                genotype.haplotype(),
                genotype.function(),
                genotype.linkedDrugs(),
                genotype.urlPrescriptionInfo(),
                genotype.panelVersion(),
                genotype.repoVersion());
    }

    private static void addCalls(@NotNull Timestamp timestamp, @NotNull InsertValuesStep12 inserter, @NotNull String sample,
            @NotNull PeachCalls calls) {
        inserter.values(timestamp,
                sample,
                calls.gene(),
                calls.positionGRCh37(),
                calls.refGRCh37(),
                calls.altGRCh37(),
                calls.positionGRCh38(),
                calls.refGRCh38(),
                calls.altGRCh38(),
                calls.rsid(),
                calls.variantAnnotation(),
                calls.filter());
    }

    void deletePeachForSample(@NotNull String sample) {
        context.delete(PEACHCALLS).where(PEACHCALLS.SAMPLEID.eq(sample)).execute();
        context.delete(PEACHGENOTYPE).where(PEACHGENOTYPE.SAMPLEID.eq(sample)).execute();
    }
}
