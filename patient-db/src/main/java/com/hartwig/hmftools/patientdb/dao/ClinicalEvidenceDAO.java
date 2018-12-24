package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CLINICALEVIDENCEITEM;

import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.actionability.EvidenceItem;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep19;

class ClinicalEvidenceDAO {

    private static final Logger LOGGER = LogManager.getLogger(ClinicalEvidenceDAO.class);

    @NotNull
    private final DSLContext context;

    ClinicalEvidenceDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeClinicalEvidence(@NotNull String sample, @NotNull List<EvidenceItem> evidenceItem) {
        deleteClinicalEvidenceForSample(sample);

        for (List<EvidenceItem> items : Iterables.partition(evidenceItem, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep19 inserter = context.insertInto(CLINICALEVIDENCEITEM,
                    CLINICALEVIDENCEITEM.SAMPLEID,
                    CLINICALEVIDENCEITEM.TYPEVARIANT,
                    CLINICALEVIDENCEITEM.GENE,
                    CLINICALEVIDENCEITEM.CHOMOSOME,
                    CLINICALEVIDENCEITEM.POSITION,
                    CLINICALEVIDENCEITEM.REF,
                    CLINICALEVIDENCEITEM.ALT,
                    CLINICALEVIDENCEITEM.CNVTYPE,
                    CLINICALEVIDENCEITEM.FUSIONFIVEGENE,
                    CLINICALEVIDENCEITEM.FUSIONTHREEGENE,
                    CLINICALEVIDENCEITEM.EVENTTYPE,
                    CLINICALEVIDENCEITEM.EVENTMATCH,
                    CLINICALEVIDENCEITEM.DRUG,
                    CLINICALEVIDENCEITEM.DRUGSTYPE,
                    CLINICALEVIDENCEITEM.RESPONSE,
                    CLINICALEVIDENCEITEM.CANCERTYPE,
                    CLINICALEVIDENCEITEM.LABEL,
                    CLINICALEVIDENCEITEM.EVIDENCELEVEL,
                    CLINICALEVIDENCEITEM.EVIDENCESOURCE);
            items.forEach(trial -> addValues(sample, trial, inserter));
            inserter.execute();
        }
    }

    private static void addValues(@NotNull String sample, @NotNull EvidenceItem evidenceItem, @NotNull InsertValuesStep19 inserter) {
        //noinspection unchecked
        inserter.values(sample,
                evidenceItem.type(),
                evidenceItem.gene(),
                evidenceItem.chromosome(),
                evidenceItem.position(),
                evidenceItem.ref(),
                evidenceItem.alt(),
                evidenceItem.cnvType(),
                evidenceItem.fusionFiveGene(),
                evidenceItem.fusionThreeGene(),
                evidenceItem.event(),
                evidenceItem.scope().readableString(),
                evidenceItem.drug(),
                evidenceItem.drugsType(),
                evidenceItem.response(),
                evidenceItem.cancerType(),
                evidenceItem.isOnLabel() ? "Tumor Type specific" : "Other tumor types specific",
                evidenceItem.level().readableString(),
                evidenceItem.source().sourceName());
    }

    void deleteClinicalEvidenceForSample(@NotNull String sample) {
        context.delete(CLINICALEVIDENCEITEM).where(CLINICALEVIDENCEITEM.SAMPLEID.eq(sample)).execute();
    }
}
