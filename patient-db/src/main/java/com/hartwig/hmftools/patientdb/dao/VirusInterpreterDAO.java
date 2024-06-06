package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.VIRUSANNOTATION;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusType;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep12;

public class VirusInterpreterDAO
{

    private final DSLContext context;

    VirusInterpreterDAO(final DSLContext context)
    {
        this.context = context;
    }

    void writeVirusInterpreter(final String sample, final List<AnnotatedVirus> virusAnnotations)
    {
        deleteVirusAnnotationForSample(sample);
        Timestamp timestamp = new Timestamp(new Date().getTime());

        for(List<AnnotatedVirus> virusAnnotation : Iterables.partition(virusAnnotations, DB_BATCH_INSERT_SIZE))
        {
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
                    VIRUSANNOTATION.LIKELIHOOD);
            virusAnnotation.forEach(x -> addVirusAnnotation(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addVirusAnnotation(final Timestamp timestamp, final InsertValuesStep12 inserter, final String sample,
            final AnnotatedVirus annotatedVirus)
    {
        VirusType interpretation = annotatedVirus.interpretation();
        inserter.values(timestamp,
                sample,
                annotatedVirus.taxid(),
                annotatedVirus.name(),
                annotatedVirus.qcStatus(),
                annotatedVirus.integrations(),
                interpretation == null ? null : interpretation.toString(),
                annotatedVirus.percentageCovered(),
                annotatedVirus.meanCoverage(),
                annotatedVirus.expectedClonalCoverage(),
                annotatedVirus.reported(),
                annotatedVirus.virusDriverLikelihoodType());
    }

    void deleteVirusAnnotationForSample(final String sample)
    {
        context.delete(VIRUSANNOTATION).where(VIRUSANNOTATION.SAMPLEID.eq(sample)).execute();
    }
}
