package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTBREAKEND;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTDISRUPTION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTFUSION;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collector;
import java.util.stream.Collectors;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableSimpleGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.SimpleGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnalysis;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Structuralvariantbreakend;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep14;
import org.jooq.InsertValuesStep4;
import org.jooq.InsertValuesStep5;
import org.jooq.Record;
import org.jooq.Record2;
import org.jooq.Result;
import org.jooq.types.UInteger;

public class StructuralVariantAnnotationDAO {

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantAnnotationDAO.class);

    @NotNull
    private final DSLContext context;

    public StructuralVariantAnnotationDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void deleteAnnotationsForSample(@NotNull String sampleId) {
        context.delete(STRUCTURALVARIANTFUSION).where(STRUCTURALVARIANTFUSION.SAMPLEID.eq(sampleId)).execute();
        context.delete(STRUCTURALVARIANTDISRUPTION).where(STRUCTURALVARIANTDISRUPTION.SAMPLEID.eq(sampleId)).execute();
        context.delete(STRUCTURALVARIANTBREAKEND).where(STRUCTURALVARIANTBREAKEND.SAMPLEID.eq(sampleId)).execute();
    }

    @NotNull
    public final List<SimpleGeneFusion> readGeneFusions(@NotNull final String sample) {
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

    private InsertValuesStep14 createBreakendInserter()
    {
        return context.insertInto(STRUCTURALVARIANTBREAKEND,
            STRUCTURALVARIANTBREAKEND.MODIFIED,
            STRUCTURALVARIANTBREAKEND.SAMPLEID,
            STRUCTURALVARIANTBREAKEND.ISSTARTEND,
            STRUCTURALVARIANTBREAKEND.STRUCTURALVARIANTID,
            STRUCTURALVARIANTBREAKEND.GENE,
            STRUCTURALVARIANTBREAKEND.GENEID,
            STRUCTURALVARIANTBREAKEND.TRANSCRIPTID,
            STRUCTURALVARIANTBREAKEND.ISCANONICALTRANSCRIPT,
            STRUCTURALVARIANTBREAKEND.STRAND,
            STRUCTURALVARIANTBREAKEND.EXONRANKUPSTREAM,
            STRUCTURALVARIANTBREAKEND.EXONPHASEUPSTREAM,
            STRUCTURALVARIANTBREAKEND.EXONRANKDOWNSTREAM,
            STRUCTURALVARIANTBREAKEND.EXONPHASEDOWNSTREAM,
            STRUCTURALVARIANTBREAKEND.EXONMAX);
    }

    @SuppressWarnings("unchecked")
    public void write(@NotNull StructuralVariantAnalysis analysis, @NotNull String sampleId) {
        final Timestamp timestamp = new Timestamp(new Date().getTime());

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
                        transcript.exonUpstream(),
                        transcript.exonUpstreamPhase(),
                        transcript.exonDownstream(),
                        transcript.exonDownstreamPhase(),
                        transcript.exonMax());

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
                    transcriptToDatabaseIdMap.get(disruption.linkedAnnotation())));
            disruptionInserter.execute();
        }
    }

}

