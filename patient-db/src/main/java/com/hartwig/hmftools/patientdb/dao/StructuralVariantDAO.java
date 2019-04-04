package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.INFERRED;
import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANT;

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

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStepN;
import org.jooq.Record;
import org.jooq.Record1;
import org.jooq.Result;

class StructuralVariantDAO {
    @NotNull
    private final DSLContext context;

    StructuralVariantDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @NotNull
    public final List<StructuralVariantData> read(@NotNull final String sample) {
        List<StructuralVariantData> structuralVariants = Lists.newArrayList();

        final Result<Record> result = context.select().from(STRUCTURALVARIANT).where(STRUCTURALVARIANT.SAMPLEID.eq(sample)).fetch();

        for (Record record : result) {

            boolean isSingleBreakend =
                    record.getValue(STRUCTURALVARIANT.ENDCHROMOSOME) == null && record.getValue(STRUCTURALVARIANT.ENDPOSITION) == null
                            && record.getValue(STRUCTURALVARIANT.ENDORIENTATION) == null;

            final String filterStr = record.getValue(STRUCTURALVARIANT.FILTER);

            // TEMP CAS: ploidy correction for NONE segment SVs
            Double ploidy = record.getValue(STRUCTURALVARIANT.PLOIDY);
            if(isSingleBreakend && ploidy == null && filterStr.equals(INFERRED))
            {
                ploidy = getValueNotNull(record.getValue(STRUCTURALVARIANT.ADJUSTEDSTARTCOPYNUMBERCHANGE));
            }

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
                    .ploidy(getValueNotNull(ploidy))
                    .type(isSingleBreakend
                            ? StructuralVariantType.SGL
                            : StructuralVariantType.fromAttribute(record.getValue(STRUCTURALVARIANT.TYPE)))
                    .homology(record.getValue(STRUCTURALVARIANT.STARTHOMOLOGYSEQUENCE))
                    .insertSequence(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCE))
                    .filter(filterStr)
                    .imprecise(byteToBoolean(record.getValue(STRUCTURALVARIANT.IMPRECISE)))
                    .qualityScore(record.getValue(STRUCTURALVARIANT.QUALSCORE))
                    .event(getValueNotNull(record.getValue(STRUCTURALVARIANT.EVENT)))
                    .startTumorVariantFragmentCount(getValueNotNull(record.getValue(STRUCTURALVARIANT.STARTTUMORVARIANTFRAGMENTCOUNT)))
                    .startTumorReferenceFragmentCount(getValueNotNull(record.getValue(STRUCTURALVARIANT.STARTTUMORREFERENCEFRAGMENTCOUNT)))
                    .startNormalVariantFragmentCount(getValueNotNull(record.getValue(STRUCTURALVARIANT.STARTNORMALVARIANTFRAGMENTCOUNT)))
                    .startNormalReferenceFragmentCount(getValueNotNull(record.getValue(STRUCTURALVARIANT.STARTNORMALREFERENCEFRAGMENTCOUNT)))
                    .endTumorVariantFragmentCount(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDTUMORVARIANTFRAGMENTCOUNT)))
                    .endTumorReferenceFragmentCount(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDTUMORREFERENCEFRAGMENTCOUNT)))
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
                    .startRefContext(getValueNotNull(record.getValue(STRUCTURALVARIANT.STARTREFCONTEXT)))
                    .endRefContext(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDREFCONTEXT)))
                    .insertSequenceAlignments(getValueNotNull(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCEALIGNMENTS)))
                    .insertSequenceRepeatClass(getValueNotNull(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCEREPEATCLASS)))
                    .insertSequenceRepeatType(getValueNotNull(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCEREPEATTYPE)))
                    .insertSequenceRepeatOrientation(getValueNotNull(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCEREPEATORIENTATION)))
                    .insertSequenceRepeatCoverage(getValueNotNull(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCEREPEATCOVERAGE)))
                    .recovered(byteToBoolean(record.getValue(STRUCTURALVARIANT.RECOVERED)))
                    .startAnchoringSupportDistance(getValueNotNull(record.getValue(STRUCTURALVARIANT.STARTANCHORINGSUPPORTDISTANCE)))
                    .endAnchoringSupportDistance(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDANCHORINGSUPPORTDISTANCE)))
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
                    .tumorVariantFragmentCount(record.getValue(STRUCTURALVARIANT.STARTTUMORVARIANTFRAGMENTCOUNT))
                    .tumorReferenceFragmentCount(record.getValue(STRUCTURALVARIANT.STARTTUMORREFERENCEFRAGMENTCOUNT))
                    .normalVariantFragmentCount(record.getValue(STRUCTURALVARIANT.STARTNORMALVARIANTFRAGMENTCOUNT))
                    .normalReferenceFragmentCount(record.getValue(STRUCTURALVARIANT.STARTNORMALREFERENCEFRAGMENTCOUNT))
                    .startOffset(record.getValue(STRUCTURALVARIANT.STARTINTERVALOFFSETSTART))
                    .endOffset(record.getValue(STRUCTURALVARIANT.STARTINTERVALOFFSETEND))
                    .inexactHomologyOffsetStart(record.getValue(STRUCTURALVARIANT.INEXACTHOMOLOGYOFFSETSTART))
                    .inexactHomologyOffsetEnd(record.getValue(STRUCTURALVARIANT.INEXACTHOMOLOGYOFFSETEND))
                    .refGenomeContext(record.getValue(STRUCTURALVARIANT.STARTREFCONTEXT))
                    .anchoringSupportDistance(record.getValue(STRUCTURALVARIANT.STARTANCHORINGSUPPORTDISTANCE))
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
                        .tumorVariantFragmentCount(record.getValue(STRUCTURALVARIANT.ENDTUMORVARIANTFRAGMENTCOUNT))
                        .tumorReferenceFragmentCount(record.getValue(STRUCTURALVARIANT.ENDTUMORREFERENCEFRAGMENTCOUNT))
                        .normalVariantFragmentCount(record.getValue(STRUCTURALVARIANT.ENDNORMALVARIANTFRAGMENTCOUNT))
                        .normalReferenceFragmentCount(record.getValue(STRUCTURALVARIANT.ENDNORMALREFERENCEFRAGMENTCOUNT))
                        .startOffset(record.getValue(STRUCTURALVARIANT.ENDINTERVALOFFSETSTART))
                        .endOffset(record.getValue(STRUCTURALVARIANT.ENDINTERVALOFFSETEND))
                        .refGenomeContext(record.getValue(STRUCTURALVARIANT.ENDREFCONTEXT))
                        .anchoringSupportDistance(record.getValue(STRUCTURALVARIANT.ENDANCHORINGSUPPORTDISTANCE))
                        .build();
            }

            EnrichedStructuralVariant variant = ImmutableEnrichedStructuralVariant.builder()
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
                    .insertSequenceAlignments(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCEALIGNMENTS))
                    .insertSequenceRepeatClass(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCEREPEATCLASS))
                    .insertSequenceRepeatType(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCEREPEATTYPE))
                    .insertSequenceRepeatOrientation(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCEREPEATORIENTATION))
                    .insertSequenceRepeatCoverage(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCEREPEATCOVERAGE))
                    .build();

            regions.add(variant);
        }
        return regions;
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
                    STRUCTURALVARIANT.STARTTUMORVARIANTFRAGMENTCOUNT,
                    STRUCTURALVARIANT.STARTTUMORREFERENCEFRAGMENTCOUNT,
                    STRUCTURALVARIANT.STARTNORMALVARIANTFRAGMENTCOUNT,
                    STRUCTURALVARIANT.STARTNORMALREFERENCEFRAGMENTCOUNT,
                    STRUCTURALVARIANT.ENDTUMORVARIANTFRAGMENTCOUNT,
                    STRUCTURALVARIANT.ENDTUMORREFERENCEFRAGMENTCOUNT,
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
                    STRUCTURALVARIANT.RECOVERYMETHOD,
                    STRUCTURALVARIANT.RECOVERYFILTER,
                    STRUCTURALVARIANT.STARTREFCONTEXT,
                    STRUCTURALVARIANT.ENDREFCONTEXT,
                    STRUCTURALVARIANT.INSERTSEQUENCEALIGNMENTS,
                    STRUCTURALVARIANT.INSERTSEQUENCEREPEATCLASS,
                    STRUCTURALVARIANT.INSERTSEQUENCEREPEATTYPE,
                    STRUCTURALVARIANT.INSERTSEQUENCEREPEATORIENTATION,
                    STRUCTURALVARIANT.INSERTSEQUENCEREPEATCOVERAGE,
                    STRUCTURALVARIANT.STARTANCHORINGSUPPORTDISTANCE,
                    STRUCTURALVARIANT.ENDANCHORINGSUPPORTDISTANCE,
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
                DatabaseUtil.checkStringLength(variant.start().homology(), STRUCTURALVARIANT.STARTHOMOLOGYSEQUENCE),
                variant.end() == null
                        ? null
                        : DatabaseUtil.checkStringLength(variant.end().homology(), STRUCTURALVARIANT.ENDHOMOLOGYSEQUENCE),
                DatabaseUtil.checkStringLength(variant.insertSequence(), STRUCTURALVARIANT.INSERTSEQUENCE),
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
                variant.start().tumorVariantFragmentCount(),
                variant.start().tumorReferenceFragmentCount(),
                variant.start().normalVariantFragmentCount(),
                variant.start().normalReferenceFragmentCount(),
                variant.end() == null ? null : variant.end().tumorVariantFragmentCount(),
                variant.end() == null ? null : variant.end().tumorReferenceFragmentCount(),
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
                variant.recoveryMethod(),
                variant.recoveryFilter(),
                variant.start().refGenomeContext(),
                variant.end() == null ? null : variant.end().refGenomeContext(),
                DatabaseUtil.checkStringLength(variant.insertSequenceAlignments(), STRUCTURALVARIANT.INSERTSEQUENCEALIGNMENTS),
                variant.insertSequenceRepeatClass(),
                variant.insertSequenceRepeatType(),
                variant.insertSequenceRepeatOrientation(),
                variant.insertSequenceRepeatCoverage(),
                variant.start().anchoringSupportDistance(),
                variant.end() == null ? 0 : variant.end().anchoringSupportDistance(),
                timestamp);
    }

    /*
    public void updateCalculatedPloidy(final List<StructuralVariantData> svDataList)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        List<StructuralvariantRecord> updates = Lists.newArrayList();

        try
        {
            for (final StructuralVariantData svData : svDataList)
            {
                StructuralvariantRecord update = new StructuralvariantRecord();
                update.setId(Integer.parseInt(svData.id()));
                update.setModified(timestamp);
                update.setPloidymin(DatabaseUtil.decimal(svData.ploidyMin()));
                update.setPloidymax(DatabaseUtil.decimal(svData.ploidyMax()));

                updates.add(update);
            }

            context.batchUpdate(updates).execute();
        }
        catch(Exception e)
        {
        }
    }
    */

    public void deleteStructuralVariantsForSample(@NotNull String sample) {
        context.delete(STRUCTURALVARIANT).where(STRUCTURALVARIANT.SAMPLEID.eq(sample)).execute();
    }

    private static double getValueNotNull(@Nullable Double value) {
        return value != null ? value : 0D;
    }

    private static int getValueNotNull(@Nullable Integer value) {
        return value != null ? value : 0;
    }

    private static byte getValueNotNull(@Nullable Byte value) {
        return value != null ? value : 0;
    }

    @NotNull
    private static String getValueNotNull(@Nullable String value) {
        return value != null ? value : Strings.EMPTY;
    }

    private static boolean byteToBoolean(@NotNull Byte b) {
        return b != 0;
    }
}