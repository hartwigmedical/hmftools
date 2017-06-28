package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep7;
import org.jooq.Record;
import org.jooq.Result;

class StructuralVariantDAO {
    @NotNull
    private final DSLContext context;

    StructuralVariantDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @NotNull
    List<StructuralVariant> read(@NotNull final String sample) {
        List<StructuralVariant> regions = Lists.newArrayList();

        Result<Record> result = context.select()
                .from(STRUCTURALVARIANT)
                .where(STRUCTURALVARIANT.SAMPLEID.eq(sample))
                .fetch();

        for (Record record : result) {

            StructuralVariant variant = ImmutableStructuralVariant.builder()
                    .startChromosome(record.getValue(STRUCTURALVARIANT.STARTCHROMOSOME))
                    .endChromosome(record.getValue(STRUCTURALVARIANT.ENDCHROMOSOME))
                    .startPosition(record.getValue(STRUCTURALVARIANT.STARTPOSITION))
                    .endPosition(record.getValue(STRUCTURALVARIANT.ENDPOSITION))
                    .type(StructuralVariantType.fromAttribute(record.getValue(STRUCTURALVARIANT.TYPE)))
                    .build();

            regions.add(variant);
        }

        return regions;
    }

    void write(@NotNull final String sample, @NotNull List<StructuralVariant> regions) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(STRUCTURALVARIANT).where(STRUCTURALVARIANT.SAMPLEID.eq(sample)).execute();

        for (List<StructuralVariant> splitRegions : Iterables.partition(regions, BATCH_INSERT_SIZE)) {
            InsertValuesStep7 inserter = context.insertInto(STRUCTURALVARIANT, STRUCTURALVARIANT.SAMPLEID,
                    STRUCTURALVARIANT.STARTCHROMOSOME, STRUCTURALVARIANT.ENDCHROMOSOME,
                    STRUCTURALVARIANT.STARTPOSITION, STRUCTURALVARIANT.ENDPOSITION, STRUCTURALVARIANT.TYPE,
                    STRUCTURALVARIANT.MODIFIED);
            splitRegions.forEach(x -> addRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private void addRecord(Timestamp timestamp, InsertValuesStep7 inserter, String sample, StructuralVariant region) {
        inserter.values(sample, region.startChromosome(), region.endChromosome(), region.startPosition(),
                region.endPosition(), region.type(), timestamp);
    }
}