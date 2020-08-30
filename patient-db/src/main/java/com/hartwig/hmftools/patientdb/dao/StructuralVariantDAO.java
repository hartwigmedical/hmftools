package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.INFERRED;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.getValueNotNull;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStepN;
import org.jooq.Record;
import org.jooq.Record1;
import org.jooq.Result;

class StructuralVariantDAO {

    private static final int MAX_LINKED_BY = 1024;

    @NotNull
    private final DSLContext context;

    StructuralVariantDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @NotNull
    public List<StructuralVariantData> read(@NotNull String sample) {
        List<StructuralVariantData> structuralVariants = Lists.newArrayList();

        Result<Record> result = context.select().from(STRUCTURALVARIANT).where(STRUCTURALVARIANT.SAMPLEID.eq(sample)).fetch();

        for (Record record : result) {
            StructuralVariantType type = StructuralVariantType.fromAttribute(record.getValue(STRUCTURALVARIANT.TYPE));

            String filterStr = record.getValue(STRUCTURALVARIANT.FILTER);

            if (type == SGL && filterStr.equals(INFERRED)) {
                type = INF;
            }

            boolean isSingleBreakend = (type == SGL) || (type == INF);

            // ploidy correction for NONE segment SVs
            Double ploidy = record.getValue(STRUCTURALVARIANT.JUNCTIONCOPYNUMBER);
            if (type == INF && ploidy == null) {
                ploidy = getValueNotNull(record.getValue(STRUCTURALVARIANT.ADJUSTEDCOPYNUMBERCHANGESTART));
            }

            structuralVariants.add(ImmutableStructuralVariantData.builder()
                    .id(record.getValue(STRUCTURALVARIANT.SVID))
                    .startChromosome(record.getValue(STRUCTURALVARIANT.STARTCHROMOSOME))
                    .endChromosome(isSingleBreakend ? "0" : record.getValue(STRUCTURALVARIANT.ENDCHROMOSOME))
                    .startPosition(record.getValue(STRUCTURALVARIANT.STARTPOSITION))
                    .endPosition(isSingleBreakend ? -1 : record.getValue(STRUCTURALVARIANT.ENDPOSITION))
                    .startOrientation(record.getValue(STRUCTURALVARIANT.STARTORIENTATION))
                    .endOrientation(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDORIENTATION)))
                    .startHomologySequence(record.getValue(STRUCTURALVARIANT.STARTHOMOLOGYSEQUENCE))
                    .endHomologySequence(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDHOMOLOGYSEQUENCE)))
                    .startAF(getValueNotNull(record.getValue(STRUCTURALVARIANT.STARTAF)))
                    .endAF(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDAF)))
                    .junctionCopyNumber(getValueNotNull(ploidy))
                    .adjustedStartAF(getValueNotNull(record.getValue(STRUCTURALVARIANT.ADJUSTEDAFSTART)))
                    .adjustedEndAF(getValueNotNull(record.getValue(STRUCTURALVARIANT.ADJUSTEDAFEND)))
                    .adjustedStartCopyNumber(getValueNotNull(record.getValue(STRUCTURALVARIANT.ADJUSTEDCOPYNUMBERSTART)))
                    .adjustedEndCopyNumber(getValueNotNull(record.getValue(STRUCTURALVARIANT.ADJUSTEDCOPYNUMBEREND)))
                    .adjustedStartCopyNumberChange(getValueNotNull(record.getValue(STRUCTURALVARIANT.ADJUSTEDCOPYNUMBERCHANGESTART)))
                    .adjustedEndCopyNumberChange(getValueNotNull(record.getValue(STRUCTURALVARIANT.ADJUSTEDCOPYNUMBERCHANGEEND)))
                    .insertSequence(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCE))
                    .type(type)
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
                    .vcfId(String.valueOf(record.getValue(STRUCTURALVARIANT.VCFID)))
                    .recovered(byteToBoolean(record.getValue(STRUCTURALVARIANT.RECOVERED)))
                    .recoveryMethod(getValueNotNull(record.getValue(STRUCTURALVARIANT.RECOVERYMETHOD)))
                    .recoveryFilter(getValueNotNull(record.getValue(STRUCTURALVARIANT.RECOVERYFILTER)))
                    .startRefContext(getValueNotNull(record.getValue(STRUCTURALVARIANT.STARTREFCONTEXT)))
                    .endRefContext(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDREFCONTEXT)))
                    .insertSequenceAlignments(getValueNotNull(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCEALIGNMENTS)))
                    .insertSequenceRepeatClass(getValueNotNull(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCEREPEATCLASS)))
                    .insertSequenceRepeatType(getValueNotNull(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCEREPEATTYPE)))
                    .insertSequenceRepeatOrientation(getValueNotNull(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCEREPEATORIENTATION)))
                    .insertSequenceRepeatCoverage(getValueNotNull(record.getValue(STRUCTURALVARIANT.INSERTSEQUENCEREPEATCOVERAGE)))
                    .startAnchoringSupportDistance(getValueNotNull(record.getValue(STRUCTURALVARIANT.STARTANCHORINGSUPPORTDISTANCE)))
                    .endAnchoringSupportDistance(getValueNotNull(record.getValue(STRUCTURALVARIANT.ENDANCHORINGSUPPORTDISTANCE)))
                    .build());
        }
        return structuralVariants;
    }

    @NotNull
    List<String> getSamplesList(@NotNull String sampleSearch) {
        Result<Record1<String>> result = sampleSearch.equals("")
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

    void write(@NotNull String sample, @NotNull List<StructuralVariantData> variants) {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        deleteStructuralVariantsForSample(sample);

        for (List<StructuralVariantData> batch : Iterables.partition(variants, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStepN inserter = context.insertInto(STRUCTURALVARIANT,
                    STRUCTURALVARIANT.SAMPLEID,
                    STRUCTURALVARIANT.SVID,
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
                    STRUCTURALVARIANT.ADJUSTEDAFSTART,
                    STRUCTURALVARIANT.ADJUSTEDCOPYNUMBERSTART,
                    STRUCTURALVARIANT.ADJUSTEDCOPYNUMBERCHANGESTART,
                    STRUCTURALVARIANT.ENDAF,
                    STRUCTURALVARIANT.ADJUSTEDAFEND,
                    STRUCTURALVARIANT.ADJUSTEDCOPYNUMBEREND,
                    STRUCTURALVARIANT.ADJUSTEDCOPYNUMBERCHANGEEND,
                    STRUCTURALVARIANT.JUNCTIONCOPYNUMBER,
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
            @NotNull StructuralVariantData variant) {
        boolean isSingle = variant.type() == SGL;

        inserter.values(sample,
                variant.id(),
                variant.startChromosome(),
                isSingle ? null : variant.endChromosome(),
                variant.startPosition(),
                isSingle ? null : variant.endPosition(),
                variant.startOrientation(),
                isSingle ? null : variant.endOrientation(),
                DatabaseUtil.checkStringLength(variant.startHomologySequence(), STRUCTURALVARIANT.STARTHOMOLOGYSEQUENCE),
                isSingle ? null : DatabaseUtil.checkStringLength(variant.endHomologySequence(), STRUCTURALVARIANT.ENDHOMOLOGYSEQUENCE),
                DatabaseUtil.checkStringLength(variant.insertSequence(), STRUCTURALVARIANT.INSERTSEQUENCE),
                variant.type(),
                DatabaseUtil.decimal(variant.startAF()),
                DatabaseUtil.decimal(variant.adjustedStartAF()),
                DatabaseUtil.decimal(variant.adjustedStartCopyNumber()),
                DatabaseUtil.decimal(variant.adjustedStartCopyNumberChange()),
                isSingle ? null : DatabaseUtil.decimal(variant.endAF()),
                isSingle ? null : DatabaseUtil.decimal(variant.adjustedEndAF()),
                isSingle ? null : DatabaseUtil.decimal(variant.adjustedEndCopyNumber()),
                isSingle ? null : DatabaseUtil.decimal(variant.adjustedEndCopyNumberChange()),
                variant.junctionCopyNumber(),
                variant.filter(),
                variant.imprecise(),
                DatabaseUtil.decimal(variant.qualityScore()),
                variant.event(),
                variant.startTumorVariantFragmentCount(),
                variant.startTumorReferenceFragmentCount(),
                variant.startNormalVariantFragmentCount(),
                variant.startNormalReferenceFragmentCount(),
                isSingle ? null : variant.endTumorVariantFragmentCount(),
                isSingle ? null : variant.endTumorReferenceFragmentCount(),
                isSingle ? null : variant.endNormalVariantFragmentCount(),
                isSingle ? null : variant.endNormalReferenceFragmentCount(),
                variant.startIntervalOffsetStart(),
                variant.startIntervalOffsetEnd(),
                isSingle ? null : variant.endIntervalOffsetStart(),
                isSingle ? null : variant.endIntervalOffsetEnd(),
                variant.inexactHomologyOffsetStart(),
                variant.inexactHomologyOffsetEnd(),
                variant.id(),
                limitSizeOfCSV(MAX_LINKED_BY, variant.startLinkedBy()),
                limitSizeOfCSV(MAX_LINKED_BY, variant.endLinkedBy()),
                variant.recovered(),
                variant.recoveryMethod(),
                variant.recoveryFilter(),
                variant.startRefContext(),
                isSingle ? null : variant.endRefContext(),
                DatabaseUtil.checkStringLength(variant.insertSequenceAlignments(), STRUCTURALVARIANT.INSERTSEQUENCEALIGNMENTS),
                variant.insertSequenceRepeatClass(),
                variant.insertSequenceRepeatType(),
                variant.insertSequenceRepeatOrientation(),
                variant.insertSequenceRepeatCoverage(),
                variant.startAnchoringSupportDistance(),
                isSingle ? 0 : variant.endAnchoringSupportDistance(),
                timestamp);
    }

    void deleteStructuralVariantsForSample(@NotNull String sample) {
        context.delete(STRUCTURALVARIANT).where(STRUCTURALVARIANT.SAMPLEID.eq(sample)).execute();
    }

    private static boolean byteToBoolean(@NotNull Byte b) {
        return b != 0;
    }

    static String limitSizeOfCSV(int maxSize, String linkedBy) {
        if (linkedBy.length() <= maxSize) {
            return linkedBy;
        }

        if (!linkedBy.contains(",")) {
            return linkedBy.substring(0, maxSize);
        }

        StringJoiner joiner = new StringJoiner(",");
        String[] csv = linkedBy.split(",");
        for (String s : csv) {
            int sizeWithNewString = joiner.length() == 0
                    ? s.length()
                    : joiner.length() + 1 + s.length();

            if (sizeWithNewString > maxSize) {
                return joiner.toString();
            }
            joiner.add(s);
        }

        return joiner.toString();
    }

}