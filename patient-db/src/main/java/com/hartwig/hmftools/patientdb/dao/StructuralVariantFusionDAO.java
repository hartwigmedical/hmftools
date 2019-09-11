package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVBREAKEND;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVFUSION;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Svbreakend;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep17;
import org.jooq.InsertValuesStep19;
import org.jooq.Record;
import org.jooq.Record2;
import org.jooq.Result;
import org.jooq.types.UInteger;

public class StructuralVariantFusionDAO
{
    @NotNull
    private final DSLContext context;

    public StructuralVariantFusionDAO(@NotNull final DSLContext context)
    {
        this.context = context;
    }

    public void deleteAnnotationsForSample(@NotNull String sampleId)
    {
        context.delete(SVFUSION).where(SVFUSION.SAMPLEID.eq(sampleId)).execute();
        context.delete(SVBREAKEND).where(SVBREAKEND.SAMPLEID.eq(sampleId)).execute();
    }

    private InsertValuesStep19 createBreakendInserter()
    {
        return context.insertInto(SVBREAKEND,
                SVBREAKEND.MODIFIED,
                SVBREAKEND.SAMPLEID,
                SVBREAKEND.SVID,
                SVBREAKEND.STARTBREAKEND,
                SVBREAKEND.GENE,
                SVBREAKEND.TRANSCRIPTID,
                SVBREAKEND.CANONICALTRANSCRIPT,
                SVBREAKEND.GENEORIENTATION,
                SVBREAKEND.DISRUPTIVE,
                SVBREAKEND.REPORTEDDISRUPTION,
                SVBREAKEND.UNDISRUPTEDCOPYNUMBER,
                SVBREAKEND.REGIONTYPE,
                SVBREAKEND.CODINGCONTEXT,
                SVBREAKEND.BIOTYPE,
                SVBREAKEND.EXONICBASEPHASE,
                SVBREAKEND.NEXTSPLICEEXONRANK,
                SVBREAKEND.NEXTSPLICEEXONPHASE,
                SVBREAKEND.NEXTSPLICEDISTANCE,
                SVBREAKEND.TOTALEXONCOUNT);
    }

    /*
    public void writeBreakendsAndFusions(@NotNull String sampleId, @NotNull List<Transcript> transcripts,
            @NotNull List<GeneFusion> fusions)
    {
        context.delete(SVFUSION).where(SVFUSION.SAMPLEID.eq(sampleId)).execute();
        context.delete(SVBREAKEND).where(SVBREAKEND.SAMPLEID.eq(sampleId)).execute();

        final Timestamp timestamp = new Timestamp(new Date().getTime());

        // a map of breakend DB Ids to transcripts for the fusion DB record foreign key to the breakend table
        final Map<Transcript, Integer> transcriptToDatabaseIdMap = Maps.newHashMap();

        InsertValuesStep19 inserter = createBreakendInserter();
        List<Transcript> transcriptsList = Lists.newArrayList();

        for (int i = 0; i < transcripts.size(); ++i)
        {
            final Transcript transcript = transcripts.get(i);
            final GeneAnnotation geneAnnotation = transcript.gene();
            boolean isUpstream = transcript.isUpstream();

            inserter.values(timestamp,
                    sampleId,
                    transcript.gene().id(),
                    geneAnnotation.isStart(),
                    geneAnnotation.GeneName,
                    transcript.StableId,
                    transcript.isCanonical(),
                    isUpstream ? "Upstream" : "Downstream",
                    transcript.isDisruptive(),
                    transcript.reportableDisruption(),
                    DatabaseUtil.decimal(transcript.undisruptedCopyNumber()),
                    transcript.regionType(),
                    transcript.codingType(),
                    transcript.bioType(),
                    transcript.exactCodingBase(),
                    transcript.nextSpliceExonRank(),
                    transcript.nextSpliceExonPhase(),
                    isUpstream ? transcript.prevSpliceAcceptorDistance() : transcript.nextSpliceAcceptorDistance(),
                    transcript.ExonMax);

            transcriptsList.add(transcript);

            // batch-insert transcripts since there can be many more than the batch size per sample
            if (transcripts.size() >= DB_BATCH_INSERT_SIZE || i == transcripts.size() - 1)
            {
                @SuppressWarnings("unchecked")
                final List<UInteger> ids = inserter.returning(SVBREAKEND.ID).fetch().getValues(0, UInteger.class);

                if (ids.size() != transcriptsList.size())
                {
                    throw new RuntimeException("Not all transcripts were inserted successfully");
                }

                for (int j = 0; j < ids.size(); j++)
                {
                    transcriptToDatabaseIdMap.put(transcriptsList.get(j), ids.get(j).intValue());
                }

                inserter = createBreakendInserter();
                transcriptsList.clear();
            }
        }

        for (List<GeneFusion> batch : Iterables.partition(fusions, DB_BATCH_INSERT_SIZE))
        {
            final InsertValuesStep17 fusionInserter = context.insertInto(SVFUSION,
                    SVFUSION.MODIFIED,
                    SVFUSION.SAMPLEID,
                    SVFUSION.FIVEPRIMEBREAKENDID,
                    SVFUSION.THREEPRIMEBREAKENDID,
                    SVFUSION.NAME,
                    SVFUSION.REPORTED,
                    SVFUSION.REPORTEDTYPE,
                    SVFUSION.PHASED,
                    SVFUSION.CHAINLENGTH,
                    SVFUSION.CHAINLINKS,
                    SVFUSION.CHAINTERMINATED,
                    SVFUSION.DOMAINSKEPT,
                    SVFUSION.DOMAINSLOST,
                    SVFUSION.SKIPPEDEXONSUP,
                    SVFUSION.SKIPPEDEXONSDOWN,
                    SVFUSION.FUSEDEXONUP,
                    SVFUSION.FUSEDEXONDOWN);

            //noinspection unchecked
            batch.forEach(fusion -> fusionInserter.values(timestamp,
                    sampleId,
                    transcriptToDatabaseIdMap.get(fusion.upstreamTrans()),
                    transcriptToDatabaseIdMap.get(fusion.downstreamTrans()),
                    fusion.name(),
                    fusion.reportable(),
                    fusion.getKnownType(),
                    fusion.phaseMatched(),
                    fusion.getChainLength(),
                    fusion.getChainLinks(),
                    fusion.isTerminated(),
                    DatabaseUtil.checkStringLength(fusion.downstreamTrans().getProteinFeaturesKept(), SVFUSION.DOMAINSKEPT),
                    DatabaseUtil.checkStringLength(fusion.downstreamTrans().getProteinFeaturesLost(), SVFUSION.DOMAINSLOST),
                    fusion.getExonsSkipped(true),
                    fusion.getExonsSkipped(false),
                    fusion.getFusedExon(true),
                    fusion.getFusedExon(false)));

            fusionInserter.execute();
        }
    }
    */

    public void writeBreakendsAndFusions(@NotNull String sampleId, @NotNull List<LinxBreakend> breakends, @NotNull List<LinxFusion> fusions)
    {
        context.delete(SVFUSION).where(SVFUSION.SAMPLEID.eq(sampleId)).execute();
        context.delete(SVBREAKEND).where(SVBREAKEND.SAMPLEID.eq(sampleId)).execute();

        final Timestamp timestamp = new Timestamp(new Date().getTime());

        // a map of breakend DB Ids to transcripts for the fusion DB record foreign key to the breakend table
        final Map<Integer, Integer> breakendIdToDbIdMap = Maps.newHashMap();

        InsertValuesStep19 inserter = createBreakendInserter();
        List<LinxBreakend> insertedBreakends = Lists.newArrayList();

        for (int i = 0; i < breakends.size(); ++i)
        {
            final LinxBreakend breakend = breakends.get(i);

            inserter.values(timestamp,
                    sampleId,
                    breakend.svId(),
                    breakend.isStart(),
                    breakend.gene(),
                    breakend.transcriptId(),
                    breakend.canonical(),
                    breakend.isUpstream() ? "Upstream" : "Downstream",
                    breakend.disruptive(),
                    breakend.reportedDisruption(),
                    DatabaseUtil.decimal(breakend.undisruptedCopyNumber()),
                    breakend.regionType(),
                    breakend.codingContext(),
                    breakend.biotype(),
                    breakend.exonBasePhase(),
                    breakend.nextSpliceExonRank(),
                    breakend.nextSpliceExonPhase(),
                    breakend.nextSpliceDistance(),
                    breakend.totalExonCount());

            insertedBreakends.add(breakend);

            // batch-insert transcripts since there can be many more than the batch size per sample
            if (insertedBreakends.size() >= DB_BATCH_INSERT_SIZE || i == breakends.size() - 1)
            {
                @SuppressWarnings("unchecked")
                final List<UInteger> ids = inserter.returning(SVBREAKEND.ID).fetch().getValues(0, UInteger.class);

                if (ids.size() != insertedBreakends.size())
                {
                    throw new RuntimeException("Not all transcripts were inserted successfully");
                }

                for (int j = 0; j < ids.size(); j++)
                {
                    breakendIdToDbIdMap.put(insertedBreakends.get(j).id(), ids.get(j).intValue());
                }

                inserter = createBreakendInserter();
                insertedBreakends.clear();
            }
        }

        for (List<LinxFusion> batch : Iterables.partition(fusions, DB_BATCH_INSERT_SIZE))
        {
            final InsertValuesStep17 fusionInserter = context.insertInto(SVFUSION,
                    SVFUSION.MODIFIED,
                    SVFUSION.SAMPLEID,
                    SVFUSION.FIVEPRIMEBREAKENDID,
                    SVFUSION.THREEPRIMEBREAKENDID,
                    SVFUSION.NAME,
                    SVFUSION.REPORTED,
                    SVFUSION.REPORTEDTYPE,
                    SVFUSION.PHASED,
                    SVFUSION.CHAINLENGTH,
                    SVFUSION.CHAINLINKS,
                    SVFUSION.CHAINTERMINATED,
                    SVFUSION.DOMAINSKEPT,
                    SVFUSION.DOMAINSLOST,
                    SVFUSION.SKIPPEDEXONSUP,
                    SVFUSION.SKIPPEDEXONSDOWN,
                    SVFUSION.FUSEDEXONUP,
                    SVFUSION.FUSEDEXONDOWN);


            for(final LinxFusion fusion : batch)
            {
                Integer fivePrimeId = breakendIdToDbIdMap.get(fusion.fivePrimeBreakendId());
                Integer threePrimeId = breakendIdToDbIdMap.get(fusion.threePrimeBreakendId());

                if(fivePrimeId == null || threePrimeId == null)
                {
                    return;
                }

                fusionInserter.values(timestamp,
                        sampleId,
                        fivePrimeId,
                        threePrimeId,
                        fusion.name(),
                        fusion.reported(),
                        fusion.reportedType(),
                        fusion.phased(),
                        fusion.chainLength(),
                        fusion.chainLinks(),
                        fusion.chainTerminated(),
                        DatabaseUtil.checkStringLength(fusion.domainsKept(), SVFUSION.DOMAINSKEPT),
                        DatabaseUtil.checkStringLength(fusion.domainsLost(), SVFUSION.DOMAINSLOST),
                        fusion.skippedExonsUp(),
                        fusion.skippedExonsDown(),
                        fusion.fusedExonUp(),
                        fusion.fusedExonDown());
            }

            fusionInserter.execute();
        }
    }

    @NotNull
    public final List<ReportableGeneFusion> readGeneFusions(@NotNull final String sample)
    {
        Set<ReportableGeneFusion> fusions = Sets.newHashSet();

        Svbreakend five = SVBREAKEND.as("five");
        Svbreakend three = SVBREAKEND.as("three");
        final Result<Record2<String, String>> resultFiveGene = context.select(five.GENE, three.GENE)
                .from(SVFUSION)
                .innerJoin(five)
                .on(five.ID.eq(SVFUSION.FIVEPRIMEBREAKENDID))
                .innerJoin(three)
                .on(three.ID.eq(SVFUSION.THREEPRIMEBREAKENDID))
                .where(SVFUSION.SAMPLEID.eq(sample))
                .fetch();

        for (Record record : resultFiveGene)
        {
            fusions.add(createFusionBuilder().geneStart(record.getValue(five.GENE)).geneEnd(record.getValue(three.GENE)).build());
        }

        return Lists.newArrayList(fusions);
    }

    @NotNull
    private static ImmutableReportableGeneFusion.Builder createFusionBuilder()
    {
        return ImmutableReportableGeneFusion.builder()
                .geneContextStart(Strings.EMPTY)
                .geneContextEnd(Strings.EMPTY)
                .geneTranscriptStart(Strings.EMPTY)
                .geneTranscriptEnd(Strings.EMPTY)
                .ploidy(1D);
    }
}

