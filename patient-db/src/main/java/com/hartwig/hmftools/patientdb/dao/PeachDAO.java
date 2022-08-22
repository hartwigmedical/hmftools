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
import org.jooq.InsertValuesStep10;
import org.jooq.InsertValuesStep18;

class PeachDAO {

    @NotNull
    private final DSLContext context;

    PeachDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writePeach(@NotNull String sample, @NotNull String isolationBarcode, @NotNull List<PeachGenotype> peachGenotypes,
            @NotNull List<PeachCalls> peachCalls) {
        deletePeachForSample(sample);

        Timestamp timestamp = new Timestamp(new Date().getTime());

        for (List<PeachGenotype> genotypes : Iterables.partition(peachGenotypes, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep10 inserter = context.insertInto(PEACHGENOTYPE,
                    PEACHGENOTYPE.MODIFIED,
                    PEACHGENOTYPE.SAMPLEID,
                    PEACHGENOTYPE.ISOLATIONBARCODE,
                    PEACHGENOTYPE.GENE,
                    PEACHGENOTYPE.HAPLOTYPE,
                    PEACHGENOTYPE.FUNCTION,
                    PEACHGENOTYPE.LINKEDDRUGS,
                    PEACHGENOTYPE.URLPRESCRIPTIONINFO,
                    PEACHGENOTYPE.PANELVERSION,
                    PEACHGENOTYPE.REPOVERSION);
            genotypes.forEach(x -> addGenotype(timestamp, inserter, sample, isolationBarcode, x));
            inserter.execute();
        }

        for (List<PeachCalls> calls : Iterables.partition(peachCalls, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep18 inserter = context.insertInto(PEACHCALLS,
                    PEACHCALLS.MODIFIED,
                    PEACHCALLS.SAMPLEID,
                    PEACHCALLS.ISOLATIONBARCODE,
                    PEACHCALLS.GENE,
                    PEACHCALLS.CHROMOSOME,
                    PEACHCALLS.POSITIONV37,
                    PEACHCALLS.POSITIONV38,
                    PEACHCALLS.REFV37,
                    PEACHCALLS.REFV38,
                    PEACHCALLS.ALLELE1,
                    PEACHCALLS.ALLELE2,
                    PEACHCALLS.RSID,
                    PEACHCALLS.VARIANTANNOTATIONV37,
                    PEACHCALLS.FILTERV37,
                    PEACHCALLS.VARIANTANNOTATIONV38,
                    PEACHCALLS.FILTERV38,
                    PEACHCALLS.PANELVERSION,
                    PEACHCALLS.REPOVERSION);
            calls.forEach(x -> addCalls(timestamp, inserter, sample, isolationBarcode, x));
            inserter.execute();
        }
    }

    private static void addGenotype(@NotNull Timestamp timestamp, @NotNull InsertValuesStep10 inserter, @NotNull String sample,
            @NotNull String isolationBarcode, @NotNull PeachGenotype genotype) {
        inserter.values(timestamp,
                sample,
                isolationBarcode,
                genotype.gene(),
                genotype.haplotype(),
                genotype.function(),
                genotype.linkedDrugs(),
                genotype.urlPrescriptionInfo(),
                genotype.panelVersion(),
                genotype.repoVersion());
    }

    private static void addCalls(@NotNull Timestamp timestamp, @NotNull InsertValuesStep18 inserter, @NotNull String sample,
            @NotNull String isolationBarcode, @NotNull PeachCalls calls) {
        inserter.values(timestamp,
                sample,
                isolationBarcode,
                calls.gene(),
                calls.chromosome(),
                calls.positionV37(),
                calls.positionV38(),
                calls.refV37(),
                calls.refV38(),
                calls.allele1(),
                calls.allele2(),
                calls.rsid(),
                calls.variantAnnotationV37(),
                calls.filterV37(),
                calls.variantAnnotationV38(),
                calls.filterV38(),
                calls.panelVersion(),
                calls.repoVersion());
    }

    void deletePeachForSample(@NotNull String sample) {
        context.delete(PEACHCALLS).where(PEACHCALLS.SAMPLEID.eq(sample)).execute();
        context.delete(PEACHGENOTYPE).where(PEACHGENOTYPE.SAMPLEID.eq(sample)).execute();
    }
}
