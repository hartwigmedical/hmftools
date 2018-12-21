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

class ClinicalTrialDAO {

    private static final Logger LOGGER = LogManager.getLogger(ClinicalTrialDAO.class);

    @NotNull
    private final DSLContext context;

    ClinicalTrialDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeClinicalTrial(@NotNull String sample, @NotNull List<ClinicalTrial> clinicalTrial) {
        LOGGER.info("writeClinicalTrial in ClinicalTrialDAO class");
        LOGGER.info(sample);
        LOGGER.info(clinicalTrial);
        deleteClinicalTrialForSample(sample);

        for (List<ClinicalTrial> trials : Iterables.partition(clinicalTrial, DB_BATCH_INSERT_SIZE)) {
            LOGGER.info("Insert kolom table names");
            InsertValuesStep18 inserter = context.insertInto(CLINICALTRIALITEM,
                    CLINICALTRIALITEM.SAMPLEID,
                    CLINICALTRIALITEM.TYPEVARIANT,
                    CLINICALTRIALITEM.GENE,
                    CLINICALTRIALITEM.CHOMOSOME,
                    CLINICALTRIALITEM.POSITION,
                    CLINICALTRIALITEM.REF,
                    CLINICALTRIALITEM.ALT,
                    CLINICALTRIALITEM.CNVTYPE,
                    CLINICALTRIALITEM.FUSIONFIVEGENE,
                    CLINICALTRIALITEM.FUSIONTHREEGENE,
                    CLINICALTRIALITEM.EVENTTYPE,
                    CLINICALTRIALITEM.EVENTMATCH,
                    CLINICALTRIALITEM.TRIAL,
                    CLINICALTRIALITEM.CANCERTYPE,
                    CLINICALTRIALITEM.LABEL,
                    CLINICALTRIALITEM.CCMOID,
                    CLINICALTRIALITEM.ICLUSIONID,
                    CLINICALTRIALITEM.EVIDENCESOURCE);
            LOGGER.info("insert values");
            trials.forEach(trial -> addValues(sample, trial, inserter));
            LOGGER.info("values inserted");
            inserter.execute();
            LOGGER.info(inserter.execute());
            LOGGER.info(inserter);
            LOGGER.info("values executed");
        }
    }

    private static void addValues(@NotNull String sample, @NotNull ClinicalTrial clinicalTrial, @NotNull InsertValuesStep18 inserter) {
        //noinspection unchecked
        inserter.values(sample,
                clinicalTrial.type(),
                clinicalTrial.gene(),
                clinicalTrial.chromosome(),
                clinicalTrial.position(),
                clinicalTrial.ref(),
                clinicalTrial.alt(),
                clinicalTrial.cnvType(),
                clinicalTrial.fusionFiveGene(),
                clinicalTrial.fusionThreeGene(),
                clinicalTrial.event(),
                clinicalTrial.scope().readableString(),
                clinicalTrial.acronym(),
                clinicalTrial.cancerType(),
                clinicalTrial.isOnLabel() ? "Tumor Type specific" : "Other tumor types specific",
                CCMOId(clinicalTrial.reference()),
                IclusionId(clinicalTrial.reference()),
                clinicalTrial.source().sourceName());
        LOGGER.info(clinicalTrial.type());
        LOGGER.info(sample);
        LOGGER.info("addValues");
    }

    void deleteClinicalTrialForSample(@NotNull String sample) {
        LOGGER.info("deleteClinicalTrialForSample");
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
