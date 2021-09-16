package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.VIRUSANNOTATION;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.virus.AnnotatedVirusV1;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep8;

public class VirusInterpreterDAO {

    @NotNull
    private final DSLContext context;

    VirusInterpreterDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeVirusInterpreter(@NotNull String sample, @NotNull List<AnnotatedVirusV1> virusAnnotations) {
        deleteVirusAnnotationForSample(sample);
        Timestamp timestamp = new Timestamp(new Date().getTime());

        for (List<AnnotatedVirusV1> virusAnnotation : Iterables.partition(virusAnnotations, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep8 inserter = context.insertInto(VIRUSANNOTATION,
                    VIRUSANNOTATION.MODIFIED,
                    VIRUSANNOTATION.SAMPLEID,
                    VIRUSANNOTATION.TAXID,
                    VIRUSANNOTATION.VIRUSNAME,
                    VIRUSANNOTATION.QCSTATUS,
                    VIRUSANNOTATION.INTEGRATIONS,
                    VIRUSANNOTATION.INTERPRETATION,
                    VIRUSANNOTATION.REPORTED);
            virusAnnotation.forEach(x -> addVirusAnnotation(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addVirusAnnotation(@NotNull Timestamp timestamp, @NotNull InsertValuesStep8 inserter, @NotNull String sample,
            @NotNull AnnotatedVirusV1 annotatedVirus) {
        inserter.values(timestamp,
                sample,
                annotatedVirus.taxid(),
                annotatedVirus.name(),
                annotatedVirus.qcStatus(),
                annotatedVirus.integrations(),
                annotatedVirus.interpretation(),
                annotatedVirus.reported());
    }

    void deleteVirusAnnotationForSample(@NotNull String sample) {
        context.delete(VIRUSANNOTATION).where(VIRUSANNOTATION.SAMPLEID.eq(sample)).execute();
    }
}
