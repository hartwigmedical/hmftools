package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CPCTECRF;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CPCTECRFDATAMODEL;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRUPECRF;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRUPECRFDATAMODEL;

import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfDatamodelField;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;

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

    void clearCpct() {
        clear(CPCTECRF.getName(), CPCTECRFDATAMODEL.getName());
    }

    void clearDrup() {
        clear(DRUPECRF.getName(), DRUPECRFDATAMODEL.getName());
    }

    private void clear(@NotNull String ecrfTable, @NotNull String ecrfDatamodelTable) {
        context.truncate(ecrfTable).execute();
        context.truncate(ecrfDatamodelTable).execute();
    }

    void writeCpctDatamodel(@NotNull Iterable<EcrfDatamodelField> datamodelFields) {
        writeDatamodel(datamodelFields, this::cpctDatamodelInserter);
    }

    void writeDrupDatamodel(@NotNull Iterable<EcrfDatamodelField> datamodelFields) {
        writeDatamodel(datamodelFields, this::drupDatamodelInserter);
    }

    private void writeDatamodel(@NotNull Iterable<EcrfDatamodelField> datamodelFields, @NotNull Supplier<InsertValuesStep4> inserter) {
        context.batch(StreamSupport.stream(datamodelFields.spliterator(), false).map(field -> {
            String codeList = StringUtils.join(field.codeList().values().toArray(), ",");
            return inserter.get().values(field.name(), field.description(), codeList, field.isRelevant() ? "TRUE" : "FALSE");
        }).collect(Collectors.toList())).execute();
    }

    void writeCpctPatient(@NotNull EcrfPatient patient, boolean sequenced) {
        writePatient(patient, sequenced, this::cpctInserter);
    }

    void writeDrupPatient(@NotNull EcrfPatient patient, boolean sequenced) {
        writePatient(patient, sequenced, this::drupInserter);
    }

    private void writePatient(@NotNull EcrfPatient patient, boolean sequenced, @NotNull Supplier<InsertValuesStep14> inserter) {
        context.batch(patient.fields()
                .stream()
                .map(field -> inserter.get()
                        .values(field.patientId(),
                                field.studyEventOID(),
                                field.studyRepeatKey(),
                                field.formOID(),
                                field.formRepeatKey(),
                                field.itemGroupOID(),
                                field.itemGroupRepeatKey(),
                                field.itemOID(),
                                field.itemValue(),
                                field.status(),
                                field.locked(),
                                sequenced ? "TRUE" : "FALSE",
                                field.name(),
                                field.isRelevant() ? "TRUE" : "FALSE"))
                .collect(Collectors.toList())).execute();
    }

    @NotNull
    private InsertValuesStep4 cpctDatamodelInserter() {
        return context.insertInto(CPCTECRFDATAMODEL,
                CPCTECRFDATAMODEL.FIELDNAME,
                CPCTECRFDATAMODEL.DESCRIPTION,
                CPCTECRFDATAMODEL.CODELIST,
                CPCTECRFDATAMODEL.RELEVANT);
    }

    @NotNull
    private InsertValuesStep4 drupDatamodelInserter() {
        return context.insertInto(DRUPECRFDATAMODEL,
                DRUPECRFDATAMODEL.FIELDNAME,
                DRUPECRFDATAMODEL.DESCRIPTION,
                DRUPECRFDATAMODEL.CODELIST,
                DRUPECRFDATAMODEL.RELEVANT);
    }

    @NotNull
    private InsertValuesStep14 cpctInserter() {
        return context.insertInto(CPCTECRF,
                CPCTECRF.PATIENTID,
                CPCTECRF.STUDYEVENT,
                CPCTECRF.STUDYEVENTKEY,
                CPCTECRF.FORM,
                CPCTECRF.FORMKEY,
                CPCTECRF.ITEMGROUP,
                CPCTECRF.ITEMGROUPKEY,
                CPCTECRF.ITEM,
                CPCTECRF.ITEMVALUE,
                CPCTECRF.STATUS,
                CPCTECRF.LOCKED,
                CPCTECRF.SEQUENCED,
                CPCTECRF.FIELDNAME,
                CPCTECRF.RELEVANT);
    }

    @NotNull
    private InsertValuesStep14 drupInserter() {
        return context.insertInto(DRUPECRF,
                DRUPECRF.PATIENTID,
                DRUPECRF.STUDYEVENT,
                DRUPECRF.STUDYEVENTKEY,
                DRUPECRF.FORM,
                DRUPECRF.FORMKEY,
                DRUPECRF.ITEMGROUP,
                DRUPECRF.ITEMGROUPKEY,
                DRUPECRF.ITEM,
                DRUPECRF.ITEMVALUE,
                DRUPECRF.STATUS,
                DRUPECRF.LOCKED,
                DRUPECRF.SEQUENCED,
                DRUPECRF.FIELDNAME,
                DRUPECRF.RELEVANT);
    }
}
