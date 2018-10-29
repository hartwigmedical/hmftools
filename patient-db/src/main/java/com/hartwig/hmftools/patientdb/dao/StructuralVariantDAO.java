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
import org.jooq.InsertValuesStepN;
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

    private Double getValueNotNull(Double value) {
        return value != null ? value : 0;
    }

    private Integer getValueNotNull(Integer value) {
        return value != null ? value : 0;
    }

    private Byte getValueNotNull(Byte value) {
        return value != null ? value : 0;
    }

    private String getValueNotNull(String value) {
        return value != null ? value : "";
    }

    @NotNull
    public final List<StructuralVariantData> read(@NotNull final String sample) {
        List<StructuralVariantData> structuralVariants = Lists.newArrayList();

        final Result<Record> result = context.select().from(STRUCTURALVARIANT).where(STRUCTURALVARIANT.SAMPLEID.eq(sample)).fetch();

        for (Record record : result) {

            boolean isSingleBreakend =
                    record.getValue(STRUCTURALVARIANT.ENDCHROMOSOME) == null && record.getValue(STRUCTURALVARIANT.ENDPOSITION) == null
                            && record.getValue(STRUCTURALVARIANT.ENDORIENTATION) == null;

            structuralVariants.add(ImmutableStructuralVariantData.builder()
                    .id(String.valueOf(record.getValue(STRUCTURALVARIANT.ID)))
                    .vcfId(String.valueOf(record.getValue(STRUCTURALVARIANT.VCFID)))
                    .startChromosome(record.getValue(STRUCTURALVARIANT.STARTCHROMOSOME))
                    .endChromosome(isSingleBreakend ? "0" : record.getValue(STRUCTURALVARIANT.ENDCHROMOSOME))
                    .startPosition(record.getValue(STRUCTURALVARIANT.STARTPOSITION))
                    .endPosition(isSingleBreakend ? -1 : record.getValue(STRUCTURALVARIANT.ENDPOSITION))
                    .startOrientation(record.getValue(STRUCTURALVARIANT.STARTORIENTATION))
                    .endOrientation(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDORIENTATION)))
                    .startAF(getValueNotNull(record.getValue(STRUCTURALVARIANT.STARTAF)))
                    .adjustedStartAF(getValueNotNull(record.getValue(STRUCTURALVARIANT.ADJUSTEDSTARTAF)))
                    .adjustedStartCopyNumber(getValueNotNull(record.getValue(STRUCTURALVARIANT.ADJUSTEDSTARTCOPYNUMBER)))
                    .adjustedStartCopyNumberChange(getValueNotNull(record.getValue(STRUCTURALVARIANT.ADJUSTEDSTARTCOPYNUMBERCHANGE)))
                    .endAF(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDAF)))
                    .adjustedEndAF(getValueNotNull(record.getValue(STRUCTURALVARIANT.ADJUSTEDENDAF)))
                    .adjustedEndCopyNumber(getValueNotNull(record.getValue(STRUCTURALVARIANT.ADJUSTEDENDCOPYNUMBER)))
                    .adjustedEndCopyNumberChange(getValueNotNull(record.getValue(STRUCTURALVARIANT.ADJUSTEDENDCOPYNUMBERCHANGE)))
                    .ploidy(getValueNotNull(record.getValue(STRUCTURALVARIANT.PLOIDY)))
                    .type(isSingleBreakend
                            ? StructuralVariantType.SGL
                            : StructuralVariantType.fromAttribute(record.getValue(STRUCTURALVARIANT.TYPE)))
                    .homology(record.getValue(STRUCTURALVARIANT.STARTHOMOLOGYSEQUENCE))
                    .insertSequence(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCE))
                    .filter(record.getValue(STRUCTURALVARIANT.FILTER))
                    .imprecise(byteToBoolean(record.getValue(STRUCTURALVARIANT.IMPRECISE)))
                    .qualityScore(record.getValue(STRUCTURALVARIANT.QUALSCORE))
                    .event(record.getValue(STRUCTURALVARIANT.EVENT))
                    .startTumourVariantFragmentCount(record.getValue(STRUCTURALVARIANT.STARTTUMOURVARIANTFRAGMENTCOUNT))
                    .startTumourReferenceFragmentCount(record.getValue(STRUCTURALVARIANT.STARTTUMOURREFERENCEFRAGMENTCOUNT))
                    .startNormalVariantFragmentCount(record.getValue(STRUCTURALVARIANT.STARTNORMALVARIANTFRAGMENTCOUNT))
                    .startNormalReferenceFragmentCount(record.getValue(STRUCTURALVARIANT.STARTNORMALREFERENCEFRAGMENTCOUNT))
                    .endTumourVariantFragmentCount(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDTUMOURVARIANTFRAGMENTCOUNT)))
                    .endTumourReferenceFragmentCount(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDTUMOURREFERENCEFRAGMENTCOUNT)))
                    .endNormalVariantFragmentCount(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDNORMALVARIANTFRAGMENTCOUNT)))
                    .endNormalReferenceFragmentCount(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDNORMALREFERENCEFRAGMENTCOUNT)))
                    .startIntervalOffsetStart(getValueNotNull(record.getValue(STRUCTURALVARIANT.STARTINTERVALOFFSETSTART)))
                    .startIntervalOffsetEnd(getValueNotNull(record.getValue(STRUCTURALVARIANT.STARTINTERVALOFFSETEND)))
                    .endIntervalOffsetStart(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDINTERVALOFFSETSTART)))
                    .endIntervalOffsetEnd(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDINTERVALOFFSETEND)))
                    .inexactHomologyOffsetStart(getValueNotNull(record.getValue(STRUCTURALVARIANT.INEXACTHOMOLOGYOFFSETSTART)))
                    .inexactHomologyOffsetEnd(getValueNotNull(record.getValue(STRUCTURALVARIANT.INEXACTHOMOLOGYOFFSETEND)))
                    .startLinkedBy(getValueNotNull(record.getValue(STRUCTURALVARIANT.STARTLINKEDBY)))
                    .endLinkedBy(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDLINKEDBY)))
                    .build());
        }
        return structuralVariants;
    }

    @NotNull
    public final List<String> getSamplesList(@NotNull final String sampleSearch) {
        final Result<Record1<String>> result = sampleSearch.equals("")
                ? context.select(STRUCTURALVARIANT.SAMPLEID)
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
                    .tumourVariantFragmentCount(record.getValue(STRUCTURALVARIANT.STARTTUMOURVARIANTFRAGMENTCOUNT))
                    .tumourReferenceFragmentCount(record.getValue(STRUCTURALVARIANT.STARTTUMOURREFERENCEFRAGMENTCOUNT))
                    .normalVariantFragmentCount(record.getValue(STRUCTURALVARIANT.STARTNORMALVARIANTFRAGMENTCOUNT))
                    .normalReferenceFragmentCount(record.getValue(STRUCTURALVARIANT.STARTNORMALREFERENCEFRAGMENTCOUNT))
                    .startOffset(record.getValue(STRUCTURALVARIANT.STARTINTERVALOFFSETSTART))
                    .endOffset(record.getValue(STRUCTURALVARIANT.STARTINTERVALOFFSETEND))
                    .inexactHomologyOffsetStart(record.getValue(STRUCTURALVARIANT.INEXACTHOMOLOGYOFFSETSTART))
                    .inexactHomologyOffsetEnd(record.getValue(STRUCTURALVARIANT.INEXACTHOMOLOGYOFFSETEND))
                    .build();

            EnrichedStructuralVariantLeg end = null;
            if (record.getValue(STRUCTURALVARIANT.ENDCHROMOSOME) != null) {
                end = ImmutableEnrichedStructuralVariantLeg.builder()
                        .chromosome(record.getValue(STRUCTURALVARIANT.ENDCHROMOSOME))
                        .position(record.getValue(STRUCTURALVARIANT.ENDPOSITION))
                        .orientation(record.getValue(STRUCTURALVARIANT.ENDORIENTATION))
                        .homology(record.getValue(STRUCTURALVARIANT.ENDHOMOLOGYSEQUENCE))
                        .alleleFrequency(record.getValue(STRUCTURALVARIANT.ENDAF))
                        .adjustedAlleleFrequency(record.getValue(STRUCTURALVARIANT.ADJUSTEDENDAF))
                        .adjustedCopyNumber(record.getValue(STRUCTURALVARIANT.ADJUSTEDENDCOPYNUMBER))
                        .adjustedCopyNumberChange(record.getValue(STRUCTURALVARIANT.ADJUSTEDENDCOPYNUMBERCHANGE))
                        .tumourVariantFragmentCount(record.getValue(STRUCTURALVARIANT.ENDTUMOURVARIANTFRAGMENTCOUNT))
                        .tumourReferenceFragmentCount(record.getValue(STRUCTURALVARIANT.ENDTUMOURREFERENCEFRAGMENTCOUNT))
                        .normalVariantFragmentCount(record.getValue(STRUCTURALVARIANT.ENDNORMALVARIANTFRAGMENTCOUNT))
                        .normalReferenceFragmentCount(record.getValue(STRUCTURALVARIANT.ENDNORMALREFERENCEFRAGMENTCOUNT))
                        .startOffset(record.getValue(STRUCTURALVARIANT.ENDINTERVALOFFSETSTART))
                        .endOffset(record.getValue(STRUCTURALVARIANT.ENDINTERVALOFFSETEND))
                        .build();
            }

            final EnrichedStructuralVariant variant = ImmutableEnrichedStructuralVariant.builder()
                    .primaryKey(record.getValue(STRUCTURALVARIANT.ID))
                    .id(record.getValue(STRUCTURALVARIANT.VCFID))
                    .start(start)
                    .end(end)
                    .insertSequence(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCE))
                    .type(StructuralVariantType.fromAttribute(record.getValue(STRUCTURALVARIANT.TYPE)))
                    .ploidy(record.getValue(STRUCTURALVARIANT.PLOIDY))
                    .filter(record.getValue(STRUCTURALVARIANT.FILTER))
                    // TODO: what's the correct approach here?
                    // jooq type conversion or just manual mapping?
                    .imprecise(byteToBoolean(record.getValue(STRUCTURALVARIANT.IMPRECISE)))
                    .recovered(byteToBoolean(record.getValue(STRUCTURALVARIANT.RECOVERED)))
                    .qualityScore(record.getValue(STRUCTURALVARIANT.QUALSCORE))
                    .event(record.getValue(STRUCTURALVARIANT.EVENT))
                    .startLinkedBy(record.getValue(STRUCTURALVARIANT.STARTLINKEDBY))
                    .endLinkedBy(record.getValue(STRUCTURALVARIANT.ENDLINKEDBY))
                    .build();

            regions.add(variant);
        }
        return regions;
    }

    private static Boolean byteToBoolean(Byte b) {
        if (b == null) {
            return null;
        }
        return b != 0;
    }

    void write(@NotNull final String sample, @NotNull final List<EnrichedStructuralVariant> variants) {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        deleteStructuralVariantsForSample(sample);

        for (List<EnrichedStructuralVariant> batch : Iterables.partition(variants, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStepN inserter = context.insertInto(STRUCTURALVARIANT,
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
                    STRUCTURALVARIANT.FILTER,
                    STRUCTURALVARIANT.IMPRECISE,
                    STRUCTURALVARIANT.QUALSCORE,
                    STRUCTURALVARIANT.EVENT,
                    STRUCTURALVARIANT.STARTTUMOURVARIANTFRAGMENTCOUNT,
                    STRUCTURALVARIANT.STARTTUMOURREFERENCEFRAGMENTCOUNT,
                    STRUCTURALVARIANT.STARTNORMALVARIANTFRAGMENTCOUNT,
                    STRUCTURALVARIANT.STARTNORMALREFERENCEFRAGMENTCOUNT,
                    STRUCTURALVARIANT.ENDTUMOURVARIANTFRAGMENTCOUNT,
                    STRUCTURALVARIANT.ENDTUMOURREFERENCEFRAGMENTCOUNT,
                    STRUCTURALVARIANT.ENDNORMALVARIANTFRAGMENTCOUNT,
                    STRUCTURALVARIANT.ENDNORMALREFERENCEFRAGMENTCOUNT,
                    STRUCTURALVARIANT.STARTINTERVALOFFSETSTART,
                    STRUCTURALVARIANT.STARTINTERVALOFFSETEND,
                    STRUCTURALVARIANT.ENDINTERVALOFFSETSTART,
                    STRUCTURALVARIANT.ENDINTERVALOFFSETEND,
                    STRUCTURALVARIANT.INEXACTHOMOLOGYOFFSETSTART,
                    STRUCTURALVARIANT.INEXACTHOMOLOGYOFFSETEND,
                    STRUCTURALVARIANT.VCFID,
                    STRUCTURALVARIANT.STARTLINKEDBY,
                    STRUCTURALVARIANT.ENDLINKEDBY,
                    STRUCTURALVARIANT.RECOVERED,
                    STRUCTURALVARIANT.MODIFIED);
            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStepN inserter, @NotNull String sample,
            @NotNull EnrichedStructuralVariant variant) {
        //noinspection unchecked
        inserter.values(sample,
                variant.start().chromosome(),
                variant.end() == null ? null : variant.end().chromosome(),
                variant.start().position(),
                variant.end() == null ? null : variant.end().position(),
                variant.start().orientation(),
                variant.end() == null ? null : variant.end().orientation(),
                variant.start().homology(),
                variant.end() == null ? null : variant.end().homology(),
                variant.insertSequence(),
                variant.type(),
                DatabaseUtil.decimal(variant.start().alleleFrequency()),
                DatabaseUtil.decimal(variant.start().adjustedAlleleFrequency()),
                DatabaseUtil.decimal(variant.start().adjustedCopyNumber()),
                DatabaseUtil.decimal(variant.start().adjustedCopyNumberChange()),
                variant.end() == null ? null : DatabaseUtil.decimal(variant.end().alleleFrequency()),
                variant.end() == null ? null : DatabaseUtil.decimal(variant.end().adjustedAlleleFrequency()),
                variant.end() == null ? null : DatabaseUtil.decimal(variant.end().adjustedCopyNumber()),
                variant.end() == null ? null : DatabaseUtil.decimal(variant.end().adjustedCopyNumberChange()),
                variant.ploidy(),
                variant.filter(),
                variant.imprecise(),
                DatabaseUtil.decimal(variant.qualityScore()),
                variant.event(),
                variant.start().tumourVariantFragmentCount(),
                variant.start().tumourReferenceFragmentCount(),
                variant.start().normalVariantFragmentCount(),
                variant.start().normalReferenceFragmentCount(),
                variant.end() == null ? null : variant.end().tumourVariantFragmentCount(),
                variant.end() == null ? null : variant.end().tumourReferenceFragmentCount(),
                variant.end() == null ? null : variant.end().normalVariantFragmentCount(),
                variant.end() == null ? null : variant.end().normalReferenceFragmentCount(),
                variant.start().startOffset(),
                variant.start().endOffset(),
                variant.end() == null ? null : variant.end().startOffset(),
                variant.end() == null ? null : variant.end().endOffset(),
                variant.start().inexactHomologyOffsetStart(),
                variant.start().inexactHomologyOffsetEnd(),
                variant.id(),
                variant.startLinkedBy(),
                variant.endLinkedBy(),
                variant.recovered(),
                timestamp);
    }

    void deleteStructuralVariantsForSample(@NotNull String sample) {
        context.delete(STRUCTURALVARIANTDISRUPTION)
                .where(STRUCTURALVARIANTDISRUPTION.BREAKENDID.in(selectBreakendsForSample(sample)))
                .execute();
        context.delete(STRUCTURALVARIANTFUSION)
                .where(STRUCTURALVARIANTFUSION.FIVEPRIMEBREAKENDID.in(selectBreakendsForSample(sample)))
                .execute();
        context.delete(STRUCTURALVARIANTBREAKEND).where(STRUCTURALVARIANTBREAKEND.ID.in(selectBreakendsForSample(sample))).execute();

        context.delete(STRUCTURALVARIANT).where(STRUCTURALVARIANT.SAMPLEID.eq(sample)).execute();
    }

    @NotNull
    private List<Record1<UInteger>> selectBreakendsForSample(@NotNull String sample) {
        return context.select(STRUCTURALVARIANTBREAKEND.ID)
                .from(STRUCTURALVARIANTBREAKEND)
                .innerJoin(STRUCTURALVARIANT)
                .on(STRUCTURALVARIANT.ID.eq(STRUCTURALVARIANTBREAKEND.STRUCTURALVARIANTID))
                .where(STRUCTURALVARIANT.SAMPLEID.eq(sample))
                .fetch();
    }
}