package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTBREAKEND;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTDISRUPTION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTFUSION;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep21;
import org.jooq.Record;
import org.jooq.Record1;
import org.jooq.Result;
import org.jooq.types.UInteger;

class StructuralVariantDAO {
    @NotNull
    private final DSLContext context;

    StructuralVariantDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public final List<StructuralVariantData> read(@NotNull final String sample) {
        List<StructuralVariantData> structuralVariants = Lists.newArrayList();

        final Result<Record> result = context.select().from(STRUCTURALVARIANT).where(STRUCTURALVARIANT.SAMPLEID.eq(sample)).fetch();

        for (Record record : result) {
            structuralVariants.add(ImmutableStructuralVariantData.builder()
                    .id(String.valueOf(record.getValue(STRUCTURALVARIANT.ID)))
                    .startChromosome(record.getValue(STRUCTURALVARIANT.STARTCHROMOSOME))
                    .endChromosome(record.getValue(STRUCTURALVARIANT.ENDCHROMOSOME))
                    .startPosition(record.getValue(STRUCTURALVARIANT.STARTPOSITION))
                    .endPosition(record.getValue(STRUCTURALVARIANT.ENDPOSITION))
                    .startOrientation(record.getValue(STRUCTURALVARIANT.STARTORIENTATION))
                    .endOrientation(record.getValue(STRUCTURALVARIANT.ENDORIENTATION))
                    .startAF(record.getValue(STRUCTURALVARIANT.STARTAF))
                    .adjustedStartAF(record.getValue(STRUCTURALVARIANT.ADJUSTEDSTARTAF))
                    .adjustedStartCopyNumber(record.getValue(STRUCTURALVARIANT.ADJUSTEDSTARTCOPYNUMBER))
                    .adjustedStartCopyNumberChange(record.getValue(STRUCTURALVARIANT.ADJUSTEDSTARTCOPYNUMBERCHANGE))
                    .endAF(record.getValue(STRUCTURALVARIANT.ENDAF))
                    .adjustedEndAF(record.getValue(STRUCTURALVARIANT.ADJUSTEDENDAF))
                    .adjustedEndCopyNumber(record.getValue(STRUCTURALVARIANT.ADJUSTEDENDCOPYNUMBER))
                    .adjustedEndCopyNumberChange(record.getValue(STRUCTURALVARIANT.ADJUSTEDENDCOPYNUMBERCHANGE))
                    .ploidy(record.getValue(STRUCTURALVARIANT.PLOIDY))
                    .type(StructuralVariantType.fromAttribute(record.getValue(STRUCTURALVARIANT.TYPE)))
                    .homology(record.getValue(STRUCTURALVARIANT.STARTHOMOLOGYSEQUENCE))
                    .build());
        }

        return structuralVariants;
    }

    public final List<String> getSamplesList(@NotNull final String sampleSearch) {

        final Result<Record1<String>> result = sampleSearch.equals("") ? context.select(STRUCTURALVARIANT.SAMPLEID)
                .from(STRUCTURALVARIANT)
                .groupBy(STRUCTURALVARIANT.SAMPLEID)
                .fetch()
                : context.select(STRUCTURALVARIANT.SAMPLEID)
                        .from(STRUCTURALVARIANT)
                        .where(STRUCTURALVARIANT.SAMPLEID.like(sampleSearch))
                        .groupBy(STRUCTURALVARIANT.SAMPLEID)
                        .fetch();

        List<String> samplesList = Lists.newArrayList();

        for (Record record : result) {
            samplesList.add(record.getValue(STRUCTURALVARIANT.SAMPLEID));
        }

        return samplesList;
    }

    @NotNull
    List<EnrichedStructuralVariant> readEnrichedData(@NotNull final String sample) {
        final List<EnrichedStructuralVariant> regions = Lists.newArrayList();

        final Result<Record> result = context.select().from(STRUCTURALVARIANT).where(STRUCTURALVARIANT.SAMPLEID.eq(sample)).fetch();

        for (Record record : result) {
            final EnrichedStructuralVariantLeg start = ImmutableEnrichedStructuralVariantLeg.builder()
                    .chromosome(record.getValue(STRUCTURALVARIANT.STARTCHROMOSOME))
                    .position(record.getValue(STRUCTURALVARIANT.STARTPOSITION))
                    .orientation(record.getValue(STRUCTURALVARIANT.STARTORIENTATION))
                    .homology(record.getValue(STRUCTURALVARIANT.STARTHOMOLOGYSEQUENCE))
                    .alleleFrequency(record.getValue(STRUCTURALVARIANT.STARTAF))
                    .adjustedAlleleFrequency(record.getValue(STRUCTURALVARIANT.ADJUSTEDSTARTAF))
                    .adjustedCopyNumber(record.getValue(STRUCTURALVARIANT.ADJUSTEDSTARTCOPYNUMBER))
                    .adjustedCopyNumberChange(record.getValue(STRUCTURALVARIANT.ADJUSTEDSTARTCOPYNUMBERCHANGE))
                    .build();

            final EnrichedStructuralVariantLeg end = ImmutableEnrichedStructuralVariantLeg.builder()
                    .chromosome(record.getValue(STRUCTURALVARIANT.ENDCHROMOSOME))
                    .position(record.getValue(STRUCTURALVARIANT.ENDPOSITION))
                    .orientation(record.getValue(STRUCTURALVARIANT.ENDORIENTATION))
                    .homology(record.getValue(STRUCTURALVARIANT.ENDHOMOLOGYSEQUENCE))
                    .alleleFrequency(record.getValue(STRUCTURALVARIANT.ENDAF))
                    .adjustedAlleleFrequency(record.getValue(STRUCTURALVARIANT.ADJUSTEDENDAF))
                    .adjustedCopyNumber(record.getValue(STRUCTURALVARIANT.ADJUSTEDENDCOPYNUMBER))
                    .adjustedCopyNumberChange(record.getValue(STRUCTURALVARIANT.ADJUSTEDENDCOPYNUMBERCHANGE))
                    .build();

            final EnrichedStructuralVariant variant = ImmutableEnrichedStructuralVariant.builder()
                    .primaryKey(record.getValue(STRUCTURALVARIANT.ID))
                    .id(record.getValue(STRUCTURALVARIANT.ID).toString())
                    .start(start)
                    .end(end)
                    .insertSequence(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCE))
                    .type(StructuralVariantType.fromAttribute(record.getValue(STRUCTURALVARIANT.TYPE)))
                    .ploidy(record.getValue(STRUCTURALVARIANT.PLOIDY))
                    .build();

            regions.add(variant);
        }

        return regions;
    }

    void write(@NotNull final String sample, @NotNull final List<EnrichedStructuralVariant> variants) {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        final Result<Record1<UInteger>> breakendsToDelete = context.select(STRUCTURALVARIANTBREAKEND.ID)
                .from(STRUCTURALVARIANTBREAKEND)
                .innerJoin(STRUCTURALVARIANT)
                .on(STRUCTURALVARIANT.ID.eq(STRUCTURALVARIANTBREAKEND.STRUCTURALVARIANTID))
                .where(STRUCTURALVARIANT.SAMPLEID.eq(sample))
                .fetch();

        // first delete annotations
        context.delete(STRUCTURALVARIANTDISRUPTION).where(STRUCTURALVARIANTDISRUPTION.BREAKENDID.in(breakendsToDelete)).execute();
        context.delete(STRUCTURALVARIANTFUSION).where(STRUCTURALVARIANTFUSION.FIVEPRIMEBREAKENDID.in(breakendsToDelete)).execute();
        context.delete(STRUCTURALVARIANTBREAKEND).where(STRUCTURALVARIANTBREAKEND.ID.in(breakendsToDelete)).execute();

        // and then the structural variants
        context.delete(STRUCTURALVARIANT).where(STRUCTURALVARIANT.SAMPLEID.eq(sample)).execute();

        for (List<EnrichedStructuralVariant> batch : Iterables.partition(variants, DB_BATCH_INSERT_SIZE)) {
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
            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep21 inserter, @NotNull String sample,
            @NotNull EnrichedStructuralVariant variant) {
        //noinspection unchecked
        inserter.values(sample,
                variant.start().chromosome(),
                variant.end().chromosome(),
                variant.start().position(),
                variant.end().position(),
                variant.start().orientation(),
                variant.end().orientation(),
                variant.start().homology(),
                variant.end().homology(),
                variant.insertSequence(),
                variant.type(),
                variant.start().alleleFrequency(),
                variant.start().adjustedAlleleFrequency(),
                variant.start().adjustedCopyNumber(),
                variant.start().adjustedCopyNumberChange(),
                variant.end().alleleFrequency(),
                variant.end().adjustedAlleleFrequency(),
                variant.end().adjustedCopyNumber(),
                variant.end().adjustedCopyNumberChange(),
                variant.ploidy(),
                timestamp);
    }
}