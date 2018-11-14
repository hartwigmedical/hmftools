package com.hartwig.hmftools.svannotation.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTBREAKEND;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTDISRUPTION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTFUSION;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep13;
import org.jooq.InsertValuesStep14;
import org.jooq.InsertValuesStep2;
import org.jooq.InsertValuesStep3;
import org.jooq.InsertValuesStep4;
import org.jooq.InsertValuesStep5;
import org.jooq.types.UInteger;

public class StructuralVariantAnnotationDAO {

    @NotNull
    private final DSLContext context;

    public StructuralVariantAnnotationDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void deleteAnnotationsForSample(final String sampleId)
    {
        context.delete(STRUCTURALVARIANTFUSION).where(STRUCTURALVARIANTFUSION.SAMPLEID.eq(sampleId)).execute();

        context.delete(STRUCTURALVARIANTDISRUPTION).where(STRUCTURALVARIANTDISRUPTION.SAMPLEID.eq(sampleId)).execute();

        context.delete(STRUCTURALVARIANTBREAKEND).where(STRUCTURALVARIANTBREAKEND.SAMPLEID.eq(sampleId)).execute();
    }

    @SuppressWarnings("unchecked")
    public void write(final StructuralVariantAnalysis analysis, final String sampleId)
    {
        final Timestamp timestamp = new Timestamp(new Date().getTime());

        final Map<Transcript, Integer> id = Maps.newHashMap();

        // load transcript annotations
        for (final StructuralVariantAnnotation annotation : analysis.annotations())
        {
            for (final GeneAnnotation geneAnnotation : annotation.annotations())
            {
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

                for (final Transcript transcript : geneAnnotation.transcripts())
                {
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

                if (ids.size() != geneAnnotation.transcripts().size())
                {
                    throw new RuntimeException("not all transcripts were inserted successfully");
                }

                for (int i = 0; i < ids.size(); i++) {
                    id.put(geneAnnotation.transcripts().get(i), ids.get(i).intValue());
                }
            }
        }

        // load fusions
        final InsertValuesStep5 fusionInserter = context.insertInto(STRUCTURALVARIANTFUSION,
                STRUCTURALVARIANTFUSION.MODIFIED,
                STRUCTURALVARIANTFUSION.SAMPLEID,
                STRUCTURALVARIANTFUSION.ISREPORTED,
                STRUCTURALVARIANTFUSION.FIVEPRIMEBREAKENDID,
                STRUCTURALVARIANTFUSION.THREEPRIMEBREAKENDID);

        for (final GeneFusion fusion : analysis.fusions())
        {
            fusionInserter.values(
                    timestamp,
                    sampleId,
                    fusion.reportable(),
                    id.get(fusion.upstreamLinkedAnnotation()),
                    id.get(fusion.downstreamLinkedAnnotation()));
        }
        fusionInserter.execute();

        // load disruptions
        final InsertValuesStep4 disruptionInserter = context.insertInto(STRUCTURALVARIANTDISRUPTION,
                STRUCTURALVARIANTFUSION.MODIFIED,
                STRUCTURALVARIANTFUSION.SAMPLEID,
                STRUCTURALVARIANTDISRUPTION.ISREPORTED,
                STRUCTURALVARIANTDISRUPTION.BREAKENDID);

        for (final GeneDisruption disruption : analysis.disruptions())
        {
            disruptionInserter.values(
                    timestamp,
                    sampleId,
                    disruption.reportable(),
                    id.get(disruption.linkedAnnotation()));
        }
        disruptionInserter.execute();
    }
}

