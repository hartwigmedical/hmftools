package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CLINICALTRIALITEM;

import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.actionability.ClinicalTrial;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep18;
import org.jooq.InsertValuesStep9;

class ClinicalTrialDAO {

    private static final Logger LOGGER = LogManager.getLogger(ClinicalTrialDAO.class);

    @NotNull
    private final DSLContext context;

    ClinicalTrialDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeClinicalTrials(@NotNull String sample, @NotNull List<ClinicalTrial> clinicalTrial) {
        deleteClinicalTrialForSample(sample);

        for (List<ClinicalTrial> trials : Iterables.partition(clinicalTrial, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep9 inserter = context.insertInto(CLINICALTRIALITEM,
                    CLINICALTRIALITEM.SAMPLEID,
                    CLINICALTRIALITEM.EVENTTYPE,
                    CLINICALTRIALITEM.EVENTMATCH,
                    CLINICALTRIALITEM.TRIAL,
                    CLINICALTRIALITEM.CANCERTYPE,
                    CLINICALTRIALITEM.LABEL,
                    CLINICALTRIALITEM.CCMOID,
                    CLINICALTRIALITEM.ICLUSIONID,
                    CLINICALTRIALITEM.EVIDENCESOURCE);
            trials.forEach(trial -> addValues(sample, trial, inserter));
            inserter.execute();
        }
    }

    private static void addValues(@NotNull String sample, @NotNull ClinicalTrial clinicalTrial, @NotNull InsertValuesStep9 inserter) {
        //noinspection unchecked
        inserter.values(sample,
                clinicalTrial.event(),
                clinicalTrial.scope().readableString(),
                clinicalTrial.acronym(),
                clinicalTrial.cancerType(),
                clinicalTrial.isOnLabel() ? "Tumor Type specific" : "Other tumor types specific",
                CCMOId(clinicalTrial.reference()),
                IclusionId(clinicalTrial.reference()),
                clinicalTrial.source().sourceName());
    }

    void deleteClinicalTrialForSample(@NotNull String sample) {
        context.delete(CLINICALTRIALITEM).where(CLINICALTRIALITEM.SAMPLEID.eq(sample)).execute();
    }

    @NotNull
    private static String CCMOId(@NotNull String reference) {
        // KODU: Expected format "EXT1 (CCMO)"
        String referenceWithoutParenthesis = reference.replace(")", "");
        String[] splitExtAndCCMO = referenceWithoutParenthesis.split("\\(");
        return splitExtAndCCMO[1];
    }

    @NotNull
    private static String IclusionId(@NotNull String reference) {
        // KODU: Expected format "EXT1 (CCMO)"
        String referenceWithoutParenthesis = reference.replace(")", "");
        String[] splitExtAndCCMO = referenceWithoutParenthesis.split("\\(");
        return splitExtAndCCMO[0];
    }
}
