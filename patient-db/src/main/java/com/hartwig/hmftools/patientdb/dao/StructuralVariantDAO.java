package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTBREAKEND;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTDISRUPTION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTFUSION;

import static org.jooq.impl.DSL.select;

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
import org.jooq.InsertValuesStep14;
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
        final List<StructuralVariant> regions = Lists.newArrayList();

        final Result<Record> result = context.select().from(STRUCTURALVARIANT).where(STRUCTURALVARIANT.SAMPLEID.eq(sample)).fetch();

        for (Record record : result) {

            final StructuralVariant variant = ImmutableStructuralVariant.builder()
                    .primaryKey(record.getValue(STRUCTURALVARIANT.ID))
                    .id(record.getValue(STRUCTURALVARIANT.ID).toString())
                    .startChromosome(record.getValue(STRUCTURALVARIANT.STARTCHROMOSOME))
                    .endChromosome(record.getValue(STRUCTURALVARIANT.ENDCHROMOSOME))
                    .startPosition(record.getValue(STRUCTURALVARIANT.STARTPOSITION))
                    .endPosition(record.getValue(STRUCTURALVARIANT.ENDPOSITION))
                    .startOrientation(record.getValue(STRUCTURALVARIANT.STARTORIENTATION))
                    .endOrientation(record.getValue(STRUCTURALVARIANT.ENDORIENTATION))
                    .startHomology(record.getValue(STRUCTURALVARIANT.STARTHOMOLOGYSEQUENCE))
                    .endHomology(record.getValue(STRUCTURALVARIANT.ENDHOMOLOGYSEQUENCE))
                    .insertSequence(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCE))
                    .type(StructuralVariantType.fromAttribute(record.getValue(STRUCTURALVARIANT.TYPE)))
                    .startAF(record.getValue(STRUCTURALVARIANT.STARTAF))
                    .endAF(record.getValue(STRUCTURALVARIANT.ENDAF))
                    .build();

            regions.add(variant);
        }

        return regions;
    }

    void write(@NotNull final String sample, @NotNull final List<StructuralVariant> regions) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(STRUCTURALVARIANTDISRUPTION)
                .where(STRUCTURALVARIANTDISRUPTION.BREAKENDID.in(select(STRUCTURALVARIANTBREAKEND.ID).from(STRUCTURALVARIANTBREAKEND)
                        .innerJoin(STRUCTURALVARIANT)
                        .on(STRUCTURALVARIANT.ID.eq(STRUCTURALVARIANTBREAKEND.STRUCTURALVARIANTID))
                        .where(STRUCTURALVARIANT.SAMPLEID.eq(sample))))
                .execute();

        context.delete(STRUCTURALVARIANTFUSION)
                .where(STRUCTURALVARIANTFUSION.FIVEPRIMEBREAKENDID.in(select(STRUCTURALVARIANTBREAKEND.ID).from(STRUCTURALVARIANTBREAKEND)
                        .innerJoin(STRUCTURALVARIANT)
                        .on(STRUCTURALVARIANT.ID.eq(STRUCTURALVARIANTBREAKEND.STRUCTURALVARIANTID))
                        .where(STRUCTURALVARIANT.SAMPLEID.eq(sample))))
                .execute();

        context.delete(STRUCTURALVARIANTBREAKEND)
                .where(STRUCTURALVARIANTBREAKEND.STRUCTURALVARIANTID.in(
                        select(STRUCTURALVARIANT.ID).from(STRUCTURALVARIANT).where(STRUCTURALVARIANT.SAMPLEID.eq(sample))))
                .execute();

        context.delete(STRUCTURALVARIANT).where(STRUCTURALVARIANT.SAMPLEID.eq(sample)).execute();

        for (List<StructuralVariant> batch : Iterables.partition(regions, BATCH_INSERT_SIZE)) {
            InsertValuesStep14 inserter =
                    context.insertInto(STRUCTURALVARIANT, STRUCTURALVARIANT.SAMPLEID, STRUCTURALVARIANT.STARTCHROMOSOME,
                            STRUCTURALVARIANT.ENDCHROMOSOME, STRUCTURALVARIANT.STARTPOSITION, STRUCTURALVARIANT.ENDPOSITION,
                            STRUCTURALVARIANT.STARTORIENTATION, STRUCTURALVARIANT.ENDORIENTATION, STRUCTURALVARIANT.STARTHOMOLOGYSEQUENCE,
                            STRUCTURALVARIANT.ENDHOMOLOGYSEQUENCE, STRUCTURALVARIANT.INSERTSEQUENCE, STRUCTURALVARIANT.TYPE,
                            STRUCTURALVARIANT.STARTAF, STRUCTURALVARIANT.ENDAF, STRUCTURALVARIANT.MODIFIED);
            batch.forEach(x -> addRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private void addRecord(Timestamp timestamp, InsertValuesStep14 inserter, String sample, StructuralVariant region) {
        inserter.values(sample, region.startChromosome(), region.endChromosome(), region.startPosition(), region.endPosition(),
                region.startOrientation(), region.endOrientation(), region.startHomology(), region.endHomology(), region.insertSequence(),
                region.type(), region.startAF(), region.endAF(), timestamp);
    }
}