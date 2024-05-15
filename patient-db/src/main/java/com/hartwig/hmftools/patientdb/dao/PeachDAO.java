package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PEACHGENOTYPE;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.peach.PeachGenotype;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep9;

class PeachDAO {

    @NotNull
    private final DSLContext context;

    PeachDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writePeach(@NotNull String sample, @NotNull List<PeachGenotype> peachGenotypes) {
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

    void deletePeachForSample(@NotNull String sample) {
        context.delete(PEACHGENOTYPE).where(PEACHGENOTYPE.SAMPLEID.eq(sample)).execute();
    }
}
