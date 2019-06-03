package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation.isUpstream;
import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTBREAKEND;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTFUSION;

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
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableSimpleGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.SimpleGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnalysis;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Structuralvariantbreakend;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep19;
import org.jooq.InsertValuesStep9;
import org.jooq.Record;
import org.jooq.Record2;
import org.jooq.Result;
import org.jooq.types.UInteger;

public class StructuralVariantFusionDAO
{

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantFusionDAO.class);

    @NotNull
    private final DSLContext context;

    public StructuralVariantFusionDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void deleteAnnotationsForSample(@NotNull String sampleId)
    {
        context.delete(STRUCTURALVARIANTFUSION).where(STRUCTURALVARIANTFUSION.SAMPLEID.eq(sampleId)).execute();
        context.delete(STRUCTURALVARIANTBREAKEND).where(STRUCTURALVARIANTBREAKEND.SAMPLEID.eq(sampleId)).execute();
    }

    private InsertValuesStep19 createBreakendInserter()
    {
        return context.insertInto(STRUCTURALVARIANTBREAKEND,
            STRUCTURALVARIANTBREAKEND.MODIFIED,
            STRUCTURALVARIANTBREAKEND.SAMPLEID,
            STRUCTURALVARIANTBREAKEND.STRUCTURALVARIANTID,
            STRUCTURALVARIANTBREAKEND.STARTBREAKEND,
            STRUCTURALVARIANTBREAKEND.GENE,
            STRUCTURALVARIANTBREAKEND.GENEID,
            STRUCTURALVARIANTBREAKEND.TRANSCRIPTID,
            STRUCTURALVARIANTBREAKEND.CANONICALTRANSCRIPT,
            STRUCTURALVARIANTBREAKEND.GENEORIENTATION,
            STRUCTURALVARIANTBREAKEND.DISRUPTIVE,
            STRUCTURALVARIANTBREAKEND.REPORTEDDISRUPTION,
            STRUCTURALVARIANTBREAKEND.REGIONTYPE,
            STRUCTURALVARIANTBREAKEND.CODINGCONTEXT,
            STRUCTURALVARIANTBREAKEND.BIOTYPE,
            STRUCTURALVARIANTBREAKEND.EXACTBASEPHASE,
            STRUCTURALVARIANTBREAKEND.NEXTSPLICEEXONRANK,
            STRUCTURALVARIANTBREAKEND.NEXTSPLICEEXONPHASE,
            STRUCTURALVARIANTBREAKEND.NEXTSPLICEDISTANCE,
            STRUCTURALVARIANTBREAKEND.TOTALEXONCOUNT);
    }

    public void writeBreakendsAndFusions(@NotNull String sampleId, @NotNull List<Transcript> transcripts, @NotNull List<GeneFusion> fusions)
    {
        context.delete(STRUCTURALVARIANTFUSION).where(STRUCTURALVARIANTFUSION.SAMPLEID.eq(sampleId)).execute();
        context.delete(STRUCTURALVARIANTBREAKEND).where(STRUCTURALVARIANTBREAKEND.SAMPLEID.eq(sampleId)).execute();

        final Timestamp timestamp = new Timestamp(new Date().getTime());

        // a map of breakend DB Ids to transcripts for the fusion DB record foreign key to the breakend table
        final Map<Transcript, Integer> transcriptToDatabaseIdMap = Maps.newHashMap();

        InsertValuesStep19 inserter = createBreakendInserter();
        List<Transcript> transcriptsList = Lists.newArrayList();

        for (int i = 0; i < transcripts.size(); ++i)
        {
            final Transcript transcript = transcripts.get(i);
            final GeneAnnotation geneAnnotation = transcript.parent();
            boolean isUpstream = isUpstream(geneAnnotation);

            inserter.values(timestamp,
                    sampleId,
                    transcript.parent().id(),
                    geneAnnotation.isStart(),
                    geneAnnotation.GeneName,
                    geneAnnotation.StableId,
                    transcript.StableId,
                    transcript.isCanonical(),
                    isUpstream ? "UPSTREAM" : "DOWNSTREAM",
                    transcript.isDisruptive(),
                    false,
                    transcript.regionType(),
                    transcript.codingType(),
                    transcript.bioType(),
                    0,
                    isUpstream ? transcript.ExonUpstream : transcript.ExonDownstream,
                    isUpstream ? transcript.ExonUpstreamPhase : transcript.ExonDownstreamPhase,
                    transcript.getDistanceUpstream(),
                    transcript.ExonMax);

            transcriptsList.add(transcript);

            // batch-insert transcripts since there can be many more than the batch size per sample
            if(transcripts.size() >= DB_BATCH_INSERT_SIZE || i == transcripts.size() - 1)
            {
                final List<UInteger> ids = inserter.returning(STRUCTURALVARIANTBREAKEND.ID).fetch().getValues(0, UInteger.class);

                if (ids.size() != transcriptsList.size())
                {
                    throw new RuntimeException("not all transcripts were inserted successfully");
                }

                for (int j = 0; j < ids.size(); j++)
                {
                    transcriptToDatabaseIdMap.put(transcriptsList.get(j), ids.get(j).intValue());
                }

                inserter = createBreakendInserter();
                transcriptsList.clear();
            }
        }

        LOGGER.debug("uploading {} fusions to DB", fusions.size());

        for (List<GeneFusion> batch : Iterables.partition(fusions, DB_BATCH_INSERT_SIZE))
        {
            final InsertValuesStep9 fusionInserter = context.insertInto(STRUCTURALVARIANTFUSION,
                    STRUCTURALVARIANTFUSION.MODIFIED,
                    STRUCTURALVARIANTFUSION.SAMPLEID,
                    STRUCTURALVARIANTFUSION.FIVEPRIMEBREAKENDID,
                    STRUCTURALVARIANTFUSION.THREEPRIMEBREAKENDID,
                    STRUCTURALVARIANTFUSION.NAME,
                    STRUCTURALVARIANTFUSION.REPORTED,
                    STRUCTURALVARIANTFUSION.REPORTEDTYPE,
                    STRUCTURALVARIANTFUSION.CHAINLENGTH,
                    STRUCTURALVARIANTFUSION.SKIPPEDEXONS);

            batch.forEach(fusion -> fusionInserter.values(timestamp,
                    sampleId,
                    transcriptToDatabaseIdMap.get(fusion.upstreamTrans()),
                    transcriptToDatabaseIdMap.get(fusion.downstreamTrans()),
                    fusion.upstreamTrans().geneName() + "_" + fusion.downstreamTrans().geneName(),
                    fusion.reportable(),
                    fusion.getKnownFusionType(),
                    0,
                    fusion.getExonsSkipped(true) + fusion.getExonsSkipped(false)));

            fusionInserter.execute();
        }
    }

    @Deprecated
    public void write(@NotNull StructuralVariantAnalysis analysis, @NotNull String sampleId)
    {
        /*

        List<Transcript> transcriptList = Lists.newArrayList();

        final Map<Transcript, Integer> transcriptToDatabaseIdMap = Maps.newHashMap();

        InsertValuesStep14 inserter = createBreakendInserter();

        List<GeneAnnotation> geneAnnotations = Lists.newArrayList();

        for (StructuralVariantAnnotation annotation : analysis.annotations())
        {
            geneAnnotations.addAll(annotation.annotations());
        }

        for (int i = 0; i < geneAnnotations.size(); ++i)
        {
            GeneAnnotation geneAnnotation = geneAnnotations.get(i);

            for (final Transcript transcript : geneAnnotation.transcripts())
            {
                inserter.values(timestamp,
                        sampleId,
                        geneAnnotation.isStart(),
                        transcript.parent().id(),
                        geneAnnotation.GeneName,
                        geneAnnotation.StableId,
                        transcript.StableId,
                        transcript.isCanonical(),
                        geneAnnotation.Strand,
                        transcript.ExonUpstream,
                        transcript.ExonUpstreamPhase,
                        transcript.ExonDownstream,
                        transcript.ExonDownstreamPhase,
                        transcript.ExonMax);

                transcriptList.add(transcript);
            }

            if(transcriptList.size() >= DB_BATCH_INSERT_SIZE || i == geneAnnotations.size() - 1)
            {
                final List<UInteger> ids = inserter.returning(STRUCTURALVARIANTBREAKEND.ID).fetch().getValues(0, UInteger.class);

                if (ids.size() != transcriptList.size())
                {
                    throw new RuntimeException("not all transcripts were inserted successfully");
                }

                for (int j = 0; j < ids.size(); j++)
                {
                    transcriptToDatabaseIdMap.put(transcriptList.get(j), ids.get(j).intValue());
                }

                inserter = createBreakendInserter();
                transcriptList.clear();
            }
        }

        LOGGER.debug("uploading {} fusions to DB", analysis.fusions().size());

        for (List<GeneFusion> batch : Iterables.partition(analysis.fusions(), DB_BATCH_INSERT_SIZE))
        {
            final InsertValuesStep5 fusionInserter = context.insertInto(STRUCTURALVARIANTFUSION,
                    STRUCTURALVARIANTFUSION.MODIFIED,
                    STRUCTURALVARIANTFUSION.SAMPLEID,
                    STRUCTURALVARIANTFUSION.ISREPORTED,
                    STRUCTURALVARIANTFUSION.FIVEPRIMEBREAKENDID,
                    STRUCTURALVARIANTFUSION.THREEPRIMEBREAKENDID);

            batch.forEach(fusion -> fusionInserter.values(timestamp,
                    sampleId,
                    fusion.reportable(),
                    transcriptToDatabaseIdMap.get(fusion.upstreamTrans()),
                    transcriptToDatabaseIdMap.get(fusion.downstreamTrans())));

            fusionInserter.execute();
        }

        LOGGER.debug("uploading {} disruptions to DB", analysis.disruptions().size());

        for (List<GeneDisruption> batch : Iterables.partition(analysis.disruptions(), DB_BATCH_INSERT_SIZE))
        {
            final InsertValuesStep4 disruptionInserter = context.insertInto(STRUCTURALVARIANTDISRUPTION,
                    STRUCTURALVARIANTFUSION.MODIFIED,
                    STRUCTURALVARIANTFUSION.SAMPLEID,
                    STRUCTURALVARIANTDISRUPTION.ISREPORTED,
                    STRUCTURALVARIANTDISRUPTION.BREAKENDID);

            batch.forEach(disruption -> disruptionInserter.values(timestamp,
                    sampleId,
                    disruption.reportable(),
                    transcriptToDatabaseIdMap.get(disruption.transcript())));
            disruptionInserter.execute();
        }
        */
    }

    @NotNull
    public final List<SimpleGeneFusion> readGeneFusions(@NotNull final String sample)
    {
        Set<SimpleGeneFusion> simpleGeneFusions = Sets.newHashSet();

        Structuralvariantbreakend five = STRUCTURALVARIANTBREAKEND.as("five");
        Structuralvariantbreakend three = STRUCTURALVARIANTBREAKEND.as("three");
        final Result<Record2<String, String>> resultFiveGene = context.select(five.GENE, three.GENE)
                .from(STRUCTURALVARIANTFUSION)
                .innerJoin(five)
                .on(five.ID.eq(STRUCTURALVARIANTFUSION.FIVEPRIMEBREAKENDID))
                .innerJoin(three)
                .on(three.ID.eq(STRUCTURALVARIANTFUSION.THREEPRIMEBREAKENDID))
                .where(STRUCTURALVARIANTFUSION.SAMPLEID.eq(sample))
                .fetch();

        for (Record record : resultFiveGene) {
            simpleGeneFusions.add(ImmutableSimpleGeneFusion.builder()
                    .fiveGene(record.getValue(five.GENE))
                    .threeGene(record.getValue(three.GENE))
                    .build());
        }

        return Lists.newArrayList(simpleGeneFusions);
    }

}

