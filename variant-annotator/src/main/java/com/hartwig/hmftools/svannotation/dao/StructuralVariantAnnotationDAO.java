package com.hartwig.hmftools.svannotation.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTBREAKEND;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTDISRUPTION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTFUSION;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;
import com.hartwig.hmftools.svannotation.annotations.GeneAnnotation;
import com.hartwig.hmftools.svannotation.annotations.GeneDisruption;
import com.hartwig.hmftools.svannotation.annotations.GeneFusion;
import com.hartwig.hmftools.svannotation.annotations.StructuralVariantAnnotation;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep13;
import org.jooq.InsertValuesStep2;
import org.jooq.InsertValuesStep3;
import org.jooq.types.UInteger;

public class StructuralVariantAnnotationDAO {

    @NotNull
    private final DSLContext context;

    public StructuralVariantAnnotationDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @SuppressWarnings("unchecked")
    public void write(final StructuralVariantAnalysis analysis) {
        final Timestamp timestamp = new Timestamp(new Date().getTime());

        final Map<Transcript, Integer> id = Maps.newHashMap();

        // NERA: load transcript annotations
        for (final StructuralVariantAnnotation annotation : analysis.annotations()) {
            for (final GeneAnnotation geneAnnotation : annotation.annotations()) {
                final InsertValuesStep13 inserter = context.insertInto(STRUCTURALVARIANTBREAKEND,
                        STRUCTURALVARIANTBREAKEND.MODIFIED,
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
                            geneAnnotation.isStart(),
                            transcript.variant().primaryKey(),
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
                    id.put(geneAnnotation.transcripts().get(i), ids.get(i).intValue());
                }
            }
        }

        // NERA: load fusions
        final InsertValuesStep3 fusionInserter = context.insertInto(STRUCTURALVARIANTFUSION,
                STRUCTURALVARIANTFUSION.ISREPORTED,
                STRUCTURALVARIANTFUSION.FIVEPRIMEBREAKENDID,
                STRUCTURALVARIANTFUSION.THREEPRIMEBREAKENDID);
        for (final GeneFusion fusion : analysis.fusions()) {
            fusionInserter.values(fusion.reportable(), id.get(fusion.upstreamLinkedAnnotation()), id.get(fusion.downstreamLinkedAnnotation()));
        }
        fusionInserter.execute();

        // NERA: load disruptions
        final InsertValuesStep2 disruptionInserter = context.insertInto(STRUCTURALVARIANTDISRUPTION,
                STRUCTURALVARIANTDISRUPTION.ISREPORTED,
                STRUCTURALVARIANTDISRUPTION.BREAKENDID);
        for (final GeneDisruption disruption : analysis.disruptions()) {
            disruptionInserter.values(disruption.reportable(), id.get(disruption.linkedAnnotation()));
        }
        disruptionInserter.execute();
    }
}

