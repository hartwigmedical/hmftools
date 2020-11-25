package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PGXCALLS;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PGXGENOTYPE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PGXVARIANT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.germline.GermlineVariant;
import com.hartwig.hmftools.common.pharmacogenetics.PGXCalls;
import com.hartwig.hmftools.common.pharmacogenetics.PGXGenotype;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Germlinevariant;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep12;
import org.jooq.InsertValuesStep9;
import org.jooq.InsertValuesStepN;

class PgxDAO {

    @NotNull
    private final DSLContext context;

    PgxDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writePgx(@NotNull String sample, @NotNull List<PGXGenotype> pgxGenotype, @NotNull List<PGXCalls> pgxCalls,
            @NotNull List<GermlineVariant> pgxVariants) {
        deletePgxForSample(sample);

        Timestamp timestamp = new Timestamp(new Date().getTime());

        for (List<PGXGenotype> genotypes : Iterables.partition(pgxGenotype, DB_BATCH_INSERT_SIZE)) {
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

        for (List<PGXCalls> calls : Iterables.partition(pgxCalls, DB_BATCH_INSERT_SIZE)) {
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

        for (List<GermlineVariant> variants : Iterables.partition(pgxVariants, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStepN inserter = context.insertInto(Tables.PGXVARIANT,
                    PGXVARIANT.SAMPLEID,
                    PGXVARIANT.CHROMOSOME,
                    PGXVARIANT.POSITION,
                    PGXVARIANT.FILTER,
                    PGXVARIANT.REFSTATUS,
                    PGXVARIANT.REPORTED,
                    PGXVARIANT.PATHOGENIC,
                    PGXVARIANT.CLINVARINFO,
                    PGXVARIANT.TYPE,
                    PGXVARIANT.REF,
                    PGXVARIANT.ALT,
                    PGXVARIANT.GENE,
                    PGXVARIANT.TRANSCRIPT,
                    PGXVARIANT.EFFECT,
                    PGXVARIANT.CODINGEFFECT,
                    PGXVARIANT.HGVSCODING,
                    PGXVARIANT.HGVSPROTEIN,
                    PGXVARIANT.MICROHOMOLOGY,
                    PGXVARIANT.REPEATSEQUENCE,
                    PGXVARIANT.REPEATCOUNT,
                    PGXVARIANT.TRINUCLEOTIDECONTEXT,
                    PGXVARIANT.ALLELEREADCOUNT,
                    PGXVARIANT.TOTALREADCOUNT,
                    PGXVARIANT.ADJUSTEDCOPYNUMBER,
                    PGXVARIANT.MINORALLELECOPYNUMBER,
                    PGXVARIANT.ADJUSTEDVAF,
                    PGXVARIANT.BIALLELIC,
                    PGXVARIANT.MODIFIED);

            variants.forEach(x -> addVariants(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addGenotype(@NotNull Timestamp timestamp, @NotNull InsertValuesStep9 inserter, @NotNull String sample,
            @NotNull PGXGenotype genotype) {
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
            @NotNull PGXCalls calls) {
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

    private static void addVariants(@NotNull Timestamp timestamp, @NotNull InsertValuesStepN inserter, @NotNull String sample,
            @NotNull GermlineVariant variant) {
        inserter.values(sample,
                variant.chromosome(),
                variant.position(),
                variant.filter(),
                variant.refStatus(),
                variant.reported(),
                variant.pathogenic(),
                variant.clinvarInfo(),
                variant.type(),
                variant.ref(),
                variant.alts(),
                variant.gene(),
                variant.transcriptId(),
                variant.effects(),
                variant.codingEffect(),
                variant.hgvsCoding(),
                variant.hgvsProtein(),
                variant.microhomology(),
                variant.repeatSequence(),
                variant.repeatCount(),
                variant.trinucleotideContext(),
                variant.alleleReadCount(),
                variant.totalReadCount(),
                DatabaseUtil.decimal(variant.adjustedCopyNumber()),
                DatabaseUtil.decimal(variant.minorAlleleJcn()),
                DatabaseUtil.decimal(variant.adjustedVaf()),
                variant.biallelic(),
                timestamp);
    }

    void deletePgxForSample(@NotNull String sample) {
        context.delete(PGXCALLS).where(PGXCALLS.SAMPLEID.eq(sample)).execute();
        context.delete(PGXGENOTYPE).where(PGXGENOTYPE.SAMPLEID.eq(sample)).execute();
        context.delete(PGXVARIANT).where(PGXVARIANT.SAMPLEID.eq(sample)).execute();
    }
}
