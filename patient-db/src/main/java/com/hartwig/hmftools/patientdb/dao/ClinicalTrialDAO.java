package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CLINICALTRIALS;

import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.actionability.ClinicalTrial;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep7;

public class ClinicalTrialDAO {

    @NotNull
    private final DSLContext context;

    ClinicalTrialDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeClinicalTrial(@NotNull String sample, @NotNull List<ClinicalTrial> clinicalTrial) {
        deleteClinicalTrialForSample(sample);

        for (List<ClinicalTrial> trials : Iterables.partition(clinicalTrial, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep7 inserter = context.insertInto(CLINICALTRIALS,
                    CLINICALTRIALS.SAMPLEID,
                    CLINICALTRIALS.EVENTTYPE,
                    CLINICALTRIALS.EVENTMATCH,
                    CLINICALTRIALS.TRIAL,
                    CLINICALTRIALS.CANCERTYPE,
                    CLINICALTRIALS.LABEL,
                    CLINICALTRIALS.EVIDENCESOURCE);
            trials.forEach(trial -> addValues(sample, trial, inserter));
            inserter.execute();

        }
    }

    private static void addValues(@NotNull String sample, @NotNull ClinicalTrial clinicalTrial, @NotNull InsertValuesStep7 inserter) {
        //noinspection unchecked
        inserter.values(sample,
                clinicalTrial.event(),
                clinicalTrial.scope().readableString(),
                clinicalTrial.acronym(),
                clinicalTrial.cancerType(),
                clinicalTrial.isOnLabel() ? "Tumor Type specific" : "Other tumor types specific",
                clinicalTrial.source().sourceName());

    }

    void deleteClinicalTrialForSample(@NotNull String sample) {
        context.delete(CLINICALTRIALS).where(CLINICALTRIALS.SAMPLEID.eq(sample)).execute();
    }
}
