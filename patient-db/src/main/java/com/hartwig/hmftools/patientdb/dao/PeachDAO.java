package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PGXCALLS;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PGXGENOTYPE;

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
        deletePgxForSample(sample);

        Timestamp timestamp = new Timestamp(new Date().getTime());
        //TODO rename tables

        for (List<PeachGenotype> genotypes : Iterables.partition(peachGenotypes, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep9 inserter = context.insertInto(PGXGENOTYPE,
                    PGXGENOTYPE.MODIFIED,
                    PGXGENOTYPE.SAMPLEID,
                    PGXGENOTYPE.GENE,
                    PGXGENOTYPE.HAPLOTYPE,
                    PGXGENOTYPE.FUNCTION,
                    PGXGENOTYPE.LINKEDDRUGS,
                    PGXGENOTYPE.URLPRESCRIPTIONINFO,
                    PGXGENOTYPE.PANELVERSION,
                    PGXGENOTYPE.REPOVERSION);
            genotypes.forEach(x -> addGenotype(timestamp, inserter, sample, x));
            inserter.execute();
        }

        for (List<PeachCalls> calls : Iterables.partition(peachCalls, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep12 inserter = context.insertInto(PGXCALLS,
                    PGXCALLS.MODIFIED,
                    PGXCALLS.SAMPLEID,
                    PGXCALLS.GENE,
                    PGXCALLS.POSITIONGRCH37,
                    PGXCALLS.REFGRCH37,
                    PGXCALLS.ALTGRCH37,
                    PGXCALLS.POSITIONGRCH38,
                    PGXCALLS.REFGRCH38,
                    PGXCALLS.ALTGRCH38,
                    PGXCALLS.RSID,
                    PGXCALLS.VARIANTANNOTATION,
                    PGXCALLS.FILTER);
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

    void deletePgxForSample(@NotNull String sample) {
        context.delete(PGXCALLS).where(PGXCALLS.SAMPLEID.eq(sample)).execute();
        context.delete(PGXGENOTYPE).where(PGXGENOTYPE.SAMPLEID.eq(sample)).execute();
    }
}
