package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRUPECRF;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRUPECRFDATAMODEL;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.ECRF;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.ECRFDATAMODEL;

import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfDatamodelField;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep14;
import org.jooq.InsertValuesStep4;
import org.jooq.tools.StringUtils;

class EcrfDAO {

    @NotNull
    private final DSLContext context;

    EcrfDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    private void clear(@NotNull final String ecrfTable, @NotNull final String ecrfDatamodelTable) {
        context.truncate(ecrfTable).execute();
        context.truncate(ecrfDatamodelTable).execute();
    }

    void clearCpct() {
        clear(ECRF.getName(), ECRFDATAMODEL.getName());
    }

    void clearDrup() {
        clear(DRUPECRF.getName(), DRUPECRFDATAMODEL.getName());
    }

    void writeCpctDatamodel(@NotNull final Iterable<EcrfDatamodelField> datamodelFields) {
        writeDatamodel(datamodelFields, this::cpctDatamodelInserter);
    }

    void writeDrupDatamodel(@NotNull final Iterable<EcrfDatamodelField> datamodelFields) {
        writeDatamodel(datamodelFields, this::drupDatamodelInserter);
    }

    void writeCpctPatient(@NotNull final EcrfPatient patient, final boolean sequenced) {
        writePatient(patient, sequenced, this::cpctInserter);
    }

    void writeDrupPatient(@NotNull final EcrfPatient patient, final boolean sequenced) {
        writePatient(patient, sequenced, this::drupInserter);
    }

    private void writeDatamodel(@NotNull final Iterable<EcrfDatamodelField> datamodelFields,
            @NotNull final Supplier<InsertValuesStep4> inserter) {
        context.batch(StreamSupport.stream(datamodelFields.spliterator(), false).map(field -> {
            final String codeList = StringUtils.join(field.codeList().values().toArray(), ",");
            return inserter.get().values(field.name(), field.description(), codeList, field.isRelevant() ? "TRUE" : "FALSE");
        }).collect(Collectors.toList())).execute();
    }

    private void writePatient(@NotNull final EcrfPatient patient, final boolean sequenced,
            @NotNull final Supplier<InsertValuesStep14> inserter) {
        context.batch(patient.fields()
                .stream()
                .map(field -> inserter.get()
                        .values(field.patientId(), field.studyEventOID(), field.studyRepeatKey(), field.formOID(), field.formRepeatKey(),
                                field.itemGroupOID(), field.itemGroupRepeatKey(), field.itemOID(), field.itemValue(), field.status(),
                                field.locked(), sequenced ? "TRUE" : "FALSE", field.name(), field.isRelevant() ? "TRUE" : "FALSE"))
                .collect(Collectors.toList())).execute();
    }

    private InsertValuesStep4 cpctDatamodelInserter() {
        return context.insertInto(ECRFDATAMODEL, ECRFDATAMODEL.FIELDNAME, ECRFDATAMODEL.DESCRIPTION, ECRFDATAMODEL.CODELIST,
                ECRFDATAMODEL.RELEVANT);
    }

    private InsertValuesStep4 drupDatamodelInserter() {
        return context.insertInto(DRUPECRFDATAMODEL, DRUPECRFDATAMODEL.FIELDNAME, DRUPECRFDATAMODEL.DESCRIPTION, DRUPECRFDATAMODEL.CODELIST,
                DRUPECRFDATAMODEL.RELEVANT);
    }

    private InsertValuesStep14 cpctInserter() {
        return context.insertInto(ECRF, ECRF.PATIENTID, ECRF.STUDYEVENT, ECRF.STUDYEVENTKEY, ECRF.FORM, ECRF.FORMKEY, ECRF.ITEMGROUP,
                ECRF.ITEMGROUPKEY, ECRF.ITEM, ECRF.ITEMVALUE, ECRF.STATUS, ECRF.LOCKED, ECRF.SEQUENCED, ECRF.FIELDNAME, ECRF.RELEVANT);
    }

    private InsertValuesStep14 drupInserter() {
        return context.insertInto(DRUPECRF, DRUPECRF.PATIENTID, DRUPECRF.STUDYEVENT, DRUPECRF.STUDYEVENTKEY, DRUPECRF.FORM,
                DRUPECRF.FORMKEY, DRUPECRF.ITEMGROUP, DRUPECRF.ITEMGROUPKEY, DRUPECRF.ITEM, DRUPECRF.ITEMVALUE, DRUPECRF.STATUS,
                DRUPECRF.LOCKED, DRUPECRF.SEQUENCED, DRUPECRF.FIELDNAME, DRUPECRF.RELEVANT);
    }
}
