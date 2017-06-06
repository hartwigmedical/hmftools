package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.COMPREHENSIVESOMATICVARIANT;

import java.sql.Timestamp;
import java.util.Collections;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
import org.jooq.Condition;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep11;
import org.jooq.Record;
import org.jooq.Result;

class ComprehensiveSomaticVariantDAO {
    private static final String PASS = "PASS";

    @NotNull
    private final DSLContext context;

    ComprehensiveSomaticVariantDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @NotNull
    List<SomaticVariant> read(@NotNull final String sample, boolean passOnly) {
        List<SomaticVariant> regions = Lists.newArrayList();

        Condition passCondition = passOnly
                ? COMPREHENSIVESOMATICVARIANT.FILTER.eq(PASS)
                : COMPREHENSIVESOMATICVARIANT.FILTER.isNotNull();

        Result<Record> result = context.select()
                .from(COMPREHENSIVESOMATICVARIANT)
                .where(COMPREHENSIVESOMATICVARIANT.SAMPLEID.eq(sample))
                .and(passCondition)
                .fetch();

        for (Record record : result) {

            final String ref = record.getValue(COMPREHENSIVESOMATICVARIANT.REF);
            final String alt = record.getValue(COMPREHENSIVESOMATICVARIANT.ALT);

            SomaticVariant variant = new SomaticVariant.Builder().chromosome(
                    record.getValue(COMPREHENSIVESOMATICVARIANT.CHROMOSOME))
                    .position(record.getValue(COMPREHENSIVESOMATICVARIANT.POSITION))
                    .ref(ref)
                    .alt(alt)
                    .alleleReadCount(record.getValue(COMPREHENSIVESOMATICVARIANT.ALLELEREADCOUNT))
                    .totalReadCount(record.getValue(COMPREHENSIVESOMATICVARIANT.TOTALREADCOUNT))
                    .type(VariantType.fromRefAlt(ref, alt))
                    .build();

            regions.add(variant);
        }

        Collections.sort(regions);
        return regions;
    }

    void write(@NotNull final String sample, @NotNull List<EnrichedSomaticVariant> regions) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(COMPREHENSIVESOMATICVARIANT).where(COMPREHENSIVESOMATICVARIANT.SAMPLEID.eq(sample)).execute();

        for (List<EnrichedSomaticVariant> splitRegions : Iterables.partition(regions, BATCH_INSERT_SIZE)) {
            InsertValuesStep11 inserter = context.insertInto(COMPREHENSIVESOMATICVARIANT,
                    COMPREHENSIVESOMATICVARIANT.SAMPLEID, COMPREHENSIVESOMATICVARIANT.CHROMOSOME,
                    COMPREHENSIVESOMATICVARIANT.POSITION, COMPREHENSIVESOMATICVARIANT.FILTER,
                    COMPREHENSIVESOMATICVARIANT.REF, COMPREHENSIVESOMATICVARIANT.ALT,
                    COMPREHENSIVESOMATICVARIANT.ALLELEREADCOUNT, COMPREHENSIVESOMATICVARIANT.TOTALREADCOUNT,
                    COMPREHENSIVESOMATICVARIANT.ADJUSTEDCOPYNUMBER, COMPREHENSIVESOMATICVARIANT.ADJUSTEDVAF,
                    COMPREHENSIVESOMATICVARIANT.MODIFIED);
            splitRegions.forEach(x -> addRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private void addRecord(Timestamp timestamp, InsertValuesStep11 inserter, String sample,
            EnrichedSomaticVariant region) {
        inserter.values(sample, region.chromosome(), region.position(), filter(region.filter()), region.ref(),
                region.alt(), region.alleleReadCount(), region.totalReadCount(), region.adjustedCopyNumber(),
                region.adjustedVAF(), timestamp);
    }

    private String filter(String filter) {
        return filter.equals(".") ? PASS : filter;
    }
}