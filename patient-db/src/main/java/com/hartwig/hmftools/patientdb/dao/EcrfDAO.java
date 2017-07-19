package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.ECRF;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.ECRFDATAMODEL;

import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfDatamodelField;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.tools.StringUtils;

class EcrfDAO {
    private static final Logger LOGGER = LogManager.getLogger(EcrfDAO.class);

    @NotNull
    private final DSLContext context;

    EcrfDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void clear() {
        context.truncate(ECRF).execute();
        context.truncate(ECRFDATAMODEL).execute();
    }

    void writeDatamodel(@NotNull final Iterable<EcrfDatamodelField> datamodelFields) {
        context.batch(StreamSupport.stream(datamodelFields.spliterator(), false).map(field -> {
            final String codeList = StringUtils.join(field.codeList().values().toArray(), ",");
            return context.insertInto(ECRFDATAMODEL, ECRFDATAMODEL.FIELDNAME, ECRFDATAMODEL.DESCRIPTION, ECRFDATAMODEL.CODELIST,
                    ECRFDATAMODEL.RELEVANT).values(field.name(), field.description(), codeList, field.isRelevant() ? "TRUE" : "FALSE");
        }).collect(Collectors.toList())).execute();
    }

    void writePatient(@NotNull final EcrfPatient patient, final boolean sequenced) {
        context.batch(patient.fields()
                .stream()
                .map(field -> context.insertInto(ECRF, ECRF.PATIENTID, ECRF.STUDYEVENT, ECRF.STUDYEVENTKEY, ECRF.FORM, ECRF.FORMKEY,
                        ECRF.ITEMGROUP, ECRF.ITEMGROUPKEY, ECRF.ITEM, ECRF.ITEMVALUE, ECRF.STATUS, ECRF.LOCKED, ECRF.SEQUENCED,
                        ECRF.FIELDNAME, ECRF.RELEVANT)
                        .values(field.patientId(), field.studyEventOID(), field.studyRepeatKey(), field.formOID(), field.formRepeatKey(),
                                field.itemGroupOID(), field.itemGroupRepeatKey(), field.itemOID(), field.itemValue(), field.status(),
                                field.locked(), sequenced ? "TRUE" : "FALSE", field.name(), field.isRelevant() ? "TRUE" : "FALSE"))
                .collect(Collectors.toList())).execute();
    }
}
