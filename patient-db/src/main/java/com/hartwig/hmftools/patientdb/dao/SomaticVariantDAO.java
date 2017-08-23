package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;

import java.sql.Timestamp;
import java.util.Collections;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
import org.jooq.Condition;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep21;
import org.jooq.Record;
import org.jooq.Result;

class SomaticVariantDAO {
    private static final String PASS = "PASS";

    @NotNull
    private final DSLContext context;

    SomaticVariantDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @NotNull
    List<EnrichedSomaticVariant> read(@NotNull final String sample, boolean passOnly) {
        final List<EnrichedSomaticVariant> regions = Lists.newArrayList();

        Condition passCondition = passOnly ? SOMATICVARIANT.FILTER.eq(PASS) : SOMATICVARIANT.FILTER.isNotNull();

        Result<Record> result = context.select().from(SOMATICVARIANT).where(SOMATICVARIANT.SAMPLEID.eq(sample)).and(passCondition).fetch();

        for (Record record : result) {

            final String ref = record.getValue(SOMATICVARIANT.REF);
            final String alt = record.getValue(SOMATICVARIANT.ALT);

            EnrichedSomaticVariant variant = ImmutableEnrichedSomaticVariant.builder()
                    .chromosome(record.getValue(SOMATICVARIANT.CHROMOSOME))
                    .position(record.getValue(SOMATICVARIANT.POSITION))
                    .ref(ref)
                    .alt(alt)
                    .gene(record.getValue(SOMATICVARIANT.GENE))
                    .cosmicId(record.getValue(SOMATICVARIANT.COSMICID))
                    .dbsnpId(record.getValue(SOMATICVARIANT.DBSNPID))
                    .effect(record.getValue(SOMATICVARIANT.EFFECT))
                    .filter(record.getValue(SOMATICVARIANT.FILTER))
                    .alleleReadCount(record.getValue(SOMATICVARIANT.ALLELEREADCOUNT))
                    .totalReadCount(record.getValue(SOMATICVARIANT.TOTALREADCOUNT))
                    .highConfidenceRegion(record.getValue(SOMATICVARIANT.HIGHCONFIDENCE, Boolean.class))
                    .adjustedVAF(record.getValue(SOMATICVARIANT.ADJUSTEDVAF))
                    .adjustedCopyNumber(record.getValue(SOMATICVARIANT.ADJUSTEDCOPYNUMBER))
                    .trinucleotideContext(record.getValue(SOMATICVARIANT.TRINUCLEOTIDECONTEXT))
                    .microhomology(record.getValue(SOMATICVARIANT.MICROHOMOLOGY))
                    .repeatSequence(record.getValue(SOMATICVARIANT.REPEATSEQUENCE))
                    .repeatCount(record.getValue(SOMATICVARIANT.REPEATCOUNT))
                    .refGenomeContext(record.getValue(SOMATICVARIANT.REFGENOMECONTEXT))
                    .type(VariantType.fromRefAlt(ref, alt))
                    .build();

            regions.add(variant);
        }

        Collections.sort(regions);
        return regions;
    }

    void write(@NotNull final String sample, @NotNull List<EnrichedSomaticVariant> regions) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(SOMATICVARIANT).where(SOMATICVARIANT.SAMPLEID.eq(sample)).execute();

        for (List<EnrichedSomaticVariant> splitRegions : Iterables.partition(regions, BATCH_INSERT_SIZE)) {
            InsertValuesStep21 inserter = context.insertInto(SOMATICVARIANT,
                    SOMATICVARIANT.SAMPLEID,
                    SOMATICVARIANT.CHROMOSOME,
                    SOMATICVARIANT.POSITION,
                    SOMATICVARIANT.FILTER,
                    SOMATICVARIANT.REF,
                    SOMATICVARIANT.ALT,
                    SOMATICVARIANT.GENE,
                    SOMATICVARIANT.COSMICID,
                    SOMATICVARIANT.DBSNPID,
                    SOMATICVARIANT.EFFECT,
                    SOMATICVARIANT.ALLELEREADCOUNT,
                    SOMATICVARIANT.TOTALREADCOUNT,
                    SOMATICVARIANT.ADJUSTEDCOPYNUMBER,
                    SOMATICVARIANT.ADJUSTEDVAF,
                    SOMATICVARIANT.HIGHCONFIDENCE,
                    SOMATICVARIANT.TRINUCLEOTIDECONTEXT,
                    SOMATICVARIANT.MICROHOMOLOGY,
                    SOMATICVARIANT.REPEATSEQUENCE,
                    SOMATICVARIANT.REPEATCOUNT,
                    SOMATICVARIANT.REFGENOMECONTEXT,
                    SOMATICVARIANT.MODIFIED);
            splitRegions.forEach(x -> addRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private void addRecord(Timestamp timestamp, InsertValuesStep21 inserter, String sample, EnrichedSomaticVariant region) {
        inserter.values(sample,
                region.chromosome(),
                region.position(),
                filter(region.filter()),
                region.ref(),
                region.alt(),
                region.gene(),
                region.cosmicId(),
                region.dbsnpId(),
                region.effect(),
                region.alleleReadCount(),
                region.totalReadCount(),
                region.adjustedCopyNumber(),
                region.adjustedVAF(),
                region.highConfidenceRegion(),
                region.trinucleotideContext(),
                region.microhomology(),
                region.repeatSequence(),
                region.repeatCount(),
                region.refGenomeContext(),
                timestamp);
    }

    private String filter(String filter) {
        return filter.equals(".") ? PASS : filter;
    }
}