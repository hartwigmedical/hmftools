package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation.isUpstream;
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
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableSimpleGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.SimpleGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Svbreakend;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep13;
import org.jooq.InsertValuesStep18;
import org.jooq.Record;
import org.jooq.Record2;
import org.jooq.Result;
import org.jooq.types.UInteger;

public class StructuralVariantFusionDAO
{
    @NotNull
    private final DSLContext context;

    public StructuralVariantFusionDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void deleteAnnotationsForSample(@NotNull String sampleId) {
        context.delete(SVFUSION).where(SVFUSION.SAMPLEID.eq(sampleId)).execute();
        context.delete(SVBREAKEND).where(SVBREAKEND.SAMPLEID.eq(sampleId)).execute();
    }

    private InsertValuesStep18 createBreakendInserter()
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
                SVBREAKEND.REGIONTYPE,
                SVBREAKEND.CODINGCONTEXT,
                SVBREAKEND.BIOTYPE,
                SVBREAKEND.EXONICBASEPHASE,
                SVBREAKEND.NEXTSPLICEEXONRANK,
                SVBREAKEND.NEXTSPLICEEXONPHASE,
                SVBREAKEND.NEXTSPLICEDISTANCE,
                SVBREAKEND.TOTALEXONCOUNT);
    }

    public void writeBreakendsAndFusions(@NotNull String sampleId, @NotNull List<Transcript> transcripts,
            @NotNull List<GeneFusion> fusions)
    {
        context.delete(SVFUSION).where(SVFUSION.SAMPLEID.eq(sampleId)).execute();
        context.delete(SVBREAKEND).where(SVBREAKEND.SAMPLEID.eq(sampleId)).execute();

        final Timestamp timestamp = new Timestamp(new Date().getTime());

        // a map of breakend DB Ids to transcripts for the fusion DB record foreign key to the breakend table
        final Map<Transcript, Integer> transcriptToDatabaseIdMap = Maps.newHashMap();

        InsertValuesStep18 inserter = createBreakendInserter();
        List<Transcript> transcriptsList = Lists.newArrayList();

        for (int i = 0; i < transcripts.size(); ++i)
        {
            final Transcript transcript = transcripts.get(i);
            final GeneAnnotation geneAnnotation = transcript.parent();
            boolean isUpstream = transcript.isUpstream();

            //noinspection unchecked
            inserter.values(timestamp,
                    sampleId,
                    transcript.parent().id(),
                    geneAnnotation.isStart(),
                    geneAnnotation.GeneName,
                    transcript.StableId,
                    transcript.isCanonical(),
                    isUpstream ? "Upstream" : "Downstream",
                    transcript.isDisruptive(),
                    transcript.reportableDisruption(),
                    transcript.regionType(),
                    transcript.codingType(),
                    transcript.bioType(),
                    transcript.exactCodingBase(),
                    isUpstream ? transcript.ExonUpstream : transcript.ExonDownstream,
                    isUpstream ? transcript.ExonUpstreamPhase : transcript.ExonDownstreamPhase,
                    isUpstream ? transcript.exonDistanceUp() : transcript.exonDistanceDown(),
                    transcript.ExonMax);

            transcriptsList.add(transcript);

            // batch-insert transcripts since there can be many more than the batch size per sample
            if (transcripts.size() >= DB_BATCH_INSERT_SIZE || i == transcripts.size() - 1)
            {

                @SuppressWarnings("unchecked")
                final List<UInteger> ids = inserter.returning(SVBREAKEND.ID).fetch().getValues(0, UInteger.class);

                if (ids.size() != transcriptsList.size()) {
                    throw new RuntimeException("Not all transcripts were inserted successfully");
                }

                for (int j = 0; j < ids.size(); j++) {
                    transcriptToDatabaseIdMap.put(transcriptsList.get(j), ids.get(j).intValue());
                }

                inserter = createBreakendInserter();
                transcriptsList.clear();
            }
        }

        for (List<GeneFusion> batch : Iterables.partition(fusions, DB_BATCH_INSERT_SIZE))
        {
            final InsertValuesStep13 fusionInserter = context.insertInto(SVFUSION,
                    SVFUSION.MODIFIED,
                    SVFUSION.SAMPLEID,
                    SVFUSION.FIVEPRIMEBREAKENDID,
                    SVFUSION.THREEPRIMEBREAKENDID,
                    SVFUSION.NAME,
                    SVFUSION.REPORTED,
                    SVFUSION.REPORTEDTYPE,
                    SVFUSION.CHAINLENGTH,
                    SVFUSION.CHAINLINKS,
                    SVFUSION.CHAINTERMINATED,
                    SVFUSION.DOMAINSKEPT,
                    SVFUSION.DOMAINSLOST,
                    SVFUSION.SKIPPEDEXONS);

            //noinspection unchecked
            batch.forEach(fusion -> fusionInserter.values(timestamp,
                    sampleId,
                    transcriptToDatabaseIdMap.get(fusion.upstreamTrans()),
                    transcriptToDatabaseIdMap.get(fusion.downstreamTrans()),
                    fusion.name(),
                    fusion.reportable(),
                    fusion.getKnownFusionType(),
                    fusion.getChainLength(),
                    fusion.getChainLinks(),
                    fusion.isTerminated(),
                    DatabaseUtil.checkStringLength(fusion.downstreamTrans().getProteinFeaturesKept(), SVFUSION.DOMAINSKEPT),
                    DatabaseUtil.checkStringLength(fusion.downstreamTrans().getProteinFeaturesLost(), SVFUSION.DOMAINSLOST),
                    fusion.getExonsSkipped(true) + fusion.getExonsSkipped(false)));

            fusionInserter.execute();
        }
    }

    @NotNull
    public final List<SimpleGeneFusion> readGeneFusions(@NotNull final String sample)
    {
        Set<SimpleGeneFusion> simpleGeneFusions = Sets.newHashSet();

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
            simpleGeneFusions.add(ImmutableSimpleGeneFusion.builder()
                    .fiveGene(record.getValue(five.GENE))
                    .threeGene(record.getValue(three.GENE))
                    .build());
        }

        return Lists.newArrayList(simpleGeneFusions);
    }
}

