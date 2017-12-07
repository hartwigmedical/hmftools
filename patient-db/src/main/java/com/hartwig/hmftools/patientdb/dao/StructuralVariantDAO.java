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
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantImpl;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantLegImpl;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep21;
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

            final StructuralVariantLeg start = ImmutableStructuralVariantLegImpl.builder()
                    .chromosome(record.getValue(STRUCTURALVARIANT.STARTCHROMOSOME))
                    .position(record.getValue(STRUCTURALVARIANT.STARTPOSITION))
                    .orientation(record.getValue(STRUCTURALVARIANT.STARTORIENTATION))
                    .homology(record.getValue(STRUCTURALVARIANT.STARTHOMOLOGYSEQUENCE))
                    .alleleFrequency(record.getValue(STRUCTURALVARIANT.STARTAF))
                    .build();

            final StructuralVariantLeg end = ImmutableStructuralVariantLegImpl.builder()
                    .chromosome(record.getValue(STRUCTURALVARIANT.ENDCHROMOSOME))
                    .position(record.getValue(STRUCTURALVARIANT.ENDPOSITION))
                    .orientation(record.getValue(STRUCTURALVARIANT.ENDORIENTATION))
                    .homology(record.getValue(STRUCTURALVARIANT.ENDHOMOLOGYSEQUENCE))
                    .alleleFrequency(record.getValue(STRUCTURALVARIANT.ENDAF))
                    .build();

            final StructuralVariant variant = ImmutableStructuralVariantImpl.builder()
                    .primaryKey(record.getValue(STRUCTURALVARIANT.ID))
                    .id(record.getValue(STRUCTURALVARIANT.ID).toString())
                    .start(start)
                    .end(end)
                    .insertSequence(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCE))
                    .type(StructuralVariantType.fromAttribute(record.getValue(STRUCTURALVARIANT.TYPE)))
                    .build();

            regions.add(variant);
        }

        return regions;
    }

    void write(@NotNull final String sample, @NotNull final List<EnrichedStructuralVariant> regions) {
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
                .where(STRUCTURALVARIANTBREAKEND.STRUCTURALVARIANTID.in(select(STRUCTURALVARIANT.ID).from(STRUCTURALVARIANT)
                        .where(STRUCTURALVARIANT.SAMPLEID.eq(sample))))
                .execute();

        context.delete(STRUCTURALVARIANT).where(STRUCTURALVARIANT.SAMPLEID.eq(sample)).execute();

        for (List<EnrichedStructuralVariant> batch : Iterables.partition(regions, BATCH_INSERT_SIZE)) {
            InsertValuesStep21 inserter = context.insertInto(STRUCTURALVARIANT,
                    STRUCTURALVARIANT.SAMPLEID,
                    STRUCTURALVARIANT.STARTCHROMOSOME,
                    STRUCTURALVARIANT.ENDCHROMOSOME,
                    STRUCTURALVARIANT.STARTPOSITION,
                    STRUCTURALVARIANT.ENDPOSITION,
                    STRUCTURALVARIANT.STARTORIENTATION,
                    STRUCTURALVARIANT.ENDORIENTATION,
                    STRUCTURALVARIANT.STARTHOMOLOGYSEQUENCE,
                    STRUCTURALVARIANT.ENDHOMOLOGYSEQUENCE,
                    STRUCTURALVARIANT.INSERTSEQUENCE,
                    STRUCTURALVARIANT.TYPE,
                    STRUCTURALVARIANT.STARTAF,
                    STRUCTURALVARIANT.ADJUSTEDSTARTAF,
                    STRUCTURALVARIANT.ADJUSTEDSTARTCOPYNUMBER,
                    STRUCTURALVARIANT.ADJUSTEDSTARTCOPYNUMBERCHANGE,
                    STRUCTURALVARIANT.ENDAF,
                    STRUCTURALVARIANT.ADJUSTEDENDAF,
                    STRUCTURALVARIANT.ADJUSTEDENDCOPYNUMBER,
                    STRUCTURALVARIANT.ADJUSTEDENDCOPYNUMBERCHANGE,
                    STRUCTURALVARIANT.PLOIDY,
                    STRUCTURALVARIANT.MODIFIED);
            batch.forEach(x -> addRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private void addRecord(Timestamp timestamp, InsertValuesStep21 inserter, String sample, EnrichedStructuralVariant region) {
        inserter.values(sample,
                region.start().chromosome(),
                region.end().chromosome(),
                region.start().position(),
                region.end().position(),
                region.start().orientation(),
                region.end().orientation(),
                region.start().homology(),
                region.end().homology(),
                region.insertSequence(),
                region.type(),
                region.start().alleleFrequency(),
                region.start().adjustedAlleleFrequency(),
                region.start().adjustedCopyNumber(),
                region.start().adjustedCopyNumberChange(),
                region.end().alleleFrequency(),
                region.end().adjustedAlleleFrequency(),
                region.end().adjustedCopyNumber(),
                region.end().adjustedCopyNumberChange(),
                region.ploidy(),
                timestamp);
    }
}