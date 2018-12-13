package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTBREAKEND;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTDISRUPTION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTFUSION;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnalysis;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep14;
import org.jooq.InsertValuesStep4;
import org.jooq.InsertValuesStep5;
import org.jooq.types.UInteger;

public class StructuralVariantAnnotationDAO {

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

    @SuppressWarnings("unchecked")
    public void write(@NotNull StructuralVariantAnalysis analysis, @NotNull String sampleId) {
        final Timestamp timestamp = new Timestamp(new Date().getTime());

        deleteAnnotationsForSample(sampleId);

        final Map<Transcript, Integer> transcriptToDatabaseIdMap = Maps.newHashMap();
        for (GeneAnnotation geneAnnotation : allAnnotations(analysis)) {
            final InsertValuesStep14 inserter = context.insertInto(STRUCTURALVARIANTBREAKEND,
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

            for (final Transcript transcript : geneAnnotation.transcripts()) {
                inserter.values(timestamp,
                        sampleId,
                        geneAnnotation.isStart(),
                        transcript.parent().id(),
                        geneAnnotation.geneName(),
                        geneAnnotation.stableId(),
                        transcript.transcriptId(),
                        transcript.isCanonical(),
                        geneAnnotation.strand(),
                        transcript.exonUpstream(),
                        transcript.exonUpstreamPhase(),
                        transcript.exonDownstream(),
                        transcript.exonDownstreamPhase(),
                        transcript.exonMax());
            }

            final List<UInteger> ids = inserter.returning(STRUCTURALVARIANTBREAKEND.ID).fetch().getValues(0, UInteger.class);

            if (ids.size() != geneAnnotation.transcripts().size()) {
                throw new RuntimeException("not all transcripts were inserted successfully");
            }

            for (int i = 0; i < ids.size(); i++) {
                transcriptToDatabaseIdMap.put(geneAnnotation.transcripts().get(i), ids.get(i).intValue());
            }
        }

        final InsertValuesStep5 fusionInserter = context.insertInto(STRUCTURALVARIANTFUSION,
                STRUCTURALVARIANTFUSION.MODIFIED,
                STRUCTURALVARIANTFUSION.SAMPLEID,
                STRUCTURALVARIANTFUSION.ISREPORTED,
                STRUCTURALVARIANTFUSION.FIVEPRIMEBREAKENDID,
                STRUCTURALVARIANTFUSION.THREEPRIMEBREAKENDID);

        for (List<GeneFusion> batch : Iterables.partition(analysis.fusions(), DB_BATCH_INSERT_SIZE)) {
            batch.forEach(fusion -> fusionInserter.values(timestamp,
                    sampleId,
                    fusion.reportable(),
                    transcriptToDatabaseIdMap.get(fusion.upstreamTrans()),
                    transcriptToDatabaseIdMap.get(fusion.downstreamTrans())));
            fusionInserter.execute();
        }

        final InsertValuesStep4 disruptionInserter = context.insertInto(STRUCTURALVARIANTDISRUPTION,
                STRUCTURALVARIANTFUSION.MODIFIED,
                STRUCTURALVARIANTFUSION.SAMPLEID,
                STRUCTURALVARIANTDISRUPTION.ISREPORTED,
                STRUCTURALVARIANTDISRUPTION.BREAKENDID);

        for (List<GeneDisruption> batch : Iterables.partition(analysis.disruptions(), DB_BATCH_INSERT_SIZE)) {
            batch.forEach(disruption -> disruptionInserter.values(timestamp,
                    sampleId,
                    disruption.reportable(),
                    transcriptToDatabaseIdMap.get(disruption.linkedAnnotation())));
            disruptionInserter.execute();
        }
    }

    @NotNull
    private static List<GeneAnnotation> allAnnotations(@NotNull StructuralVariantAnalysis analysis) {
        List<GeneAnnotation> geneAnnotations = Lists.newArrayList();
        for (StructuralVariantAnnotation annotation : analysis.annotations()) {
            geneAnnotations.addAll(annotation.annotations());
        }
        return geneAnnotations;
    }
}

