package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.ECRF;

import java.util.stream.Collectors;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

class EcrfDAO {
    private static final Logger LOGGER = LogManager.getLogger(EcrfDAO.class);

    @NotNull
    private final DSLContext context;

    EcrfDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void clear() {
        context.execute("SET FOREIGN_KEY_CHECKS = 0;");
        context.truncate(ECRF).execute();
        context.execute("SET FOREIGN_KEY_CHECKS = 1;");
    }

    void writePatient(@NotNull final EcrfPatient patient, final boolean sequenced) {
        LOGGER.info("writing patient: " + patient.patientId());
        context.batch(patient.fields()
                .stream()
                .map(field -> context.insertInto(ECRF, ECRF.PATIENTID, ECRF.STUDYEVENT, ECRF.STUDYEVENTKEY, ECRF.FORM, ECRF.FORMKEY,
                        ECRF.ITEMGROUP, ECRF.ITEMGROUPKEY, ECRF.ITEM, ECRF.ITEMVALUE, ECRF.STATUS, ECRF.LOCKED, ECRF.SEQUENCED)
                        .values(field.patientId(), field.studyEventOID(), field.studyRepeatKey(), field.formOID(), field.formRepeatKey(),
                                field.itemGroupOID(), field.itemGroupRepeatKey(), field.itemOID(), field.itemValue(), field.status(),
                                field.locked(), sequenced ? "TRUE" : "FALSE"))
                .collect(Collectors.toList())).execute();
        LOGGER.info("done writing patient: " + patient.patientId());
    }
}
