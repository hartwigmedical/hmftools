package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.VIRUSANNOTATION;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep11;
import org.jooq.InsertValuesStep12;

public class VirusInterpreterDAO {

    @NotNull
    private final DSLContext context;

    VirusInterpreterDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeVirusInterpreter(@NotNull String sample, @NotNull List<AnnotatedVirus> virusAnnotations) {
        deleteVirusAnnotationForSample(sample);
        Timestamp timestamp = new Timestamp(new Date().getTime());

        for (List<AnnotatedVirus> virusAnnotation : Iterables.partition(virusAnnotations, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep12 inserter = context.insertInto(VIRUSANNOTATION,
                    VIRUSANNOTATION.MODIFIED,
                    VIRUSANNOTATION.SAMPLEID,
                    VIRUSANNOTATION.TAXID,
                    VIRUSANNOTATION.VIRUSNAME,
                    VIRUSANNOTATION.QCSTATUS,
                    VIRUSANNOTATION.INTEGRATIONS,
                    VIRUSANNOTATION.INTERPRETATION,
                    VIRUSANNOTATION.PERCENTAGECOVERED,
                    VIRUSANNOTATION.MEANCOVERAGE,
                    VIRUSANNOTATION.EXPECTEDCLONALCOVERAGE,
                    VIRUSANNOTATION.REPORTED,
                    VIRUSANNOTATION.ISHIGHRISK);
            virusAnnotation.forEach(x -> addVirusAnnotation(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addVirusAnnotation(@NotNull Timestamp timestamp, @NotNull InsertValuesStep12 inserter, @NotNull String sample,
            @NotNull AnnotatedVirus annotatedVirus) {
        inserter.values(timestamp,
                sample,
                annotatedVirus.taxid(),
                annotatedVirus.name(),
                annotatedVirus.qcStatus(),
                annotatedVirus.integrations(),
                annotatedVirus.interpretation(),
                annotatedVirus.percentageCovered(),
                annotatedVirus.meanCoverage(),
                annotatedVirus.expectedClonalCoverage(),
                annotatedVirus.reported(),
                annotatedVirus.isHighRisk());
    }

    void deleteVirusAnnotationForSample(@NotNull String sample) {
        context.delete(VIRUSANNOTATION).where(VIRUSANNOTATION.SAMPLEID.eq(sample)).execute();
    }
}
