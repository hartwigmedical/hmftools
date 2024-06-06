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

class EcrfDAO
{
    private final DSLContext context;

    EcrfDAO(final DSLContext context)
    {
        this.context = context;
    }

    void clearCpct()
    {
        clear(CPCTECRF.getName(), CPCTECRFDATAMODEL.getName());
    }

    void clearDrup()
    {
        clear(DRUPECRF.getName(), DRUPECRFDATAMODEL.getName());
    }

    private void clear(final String ecrfTable, final String ecrfDatamodelTable)
    {
        context.truncate(ecrfTable).execute();
        context.truncate(ecrfDatamodelTable).execute();
    }

    void writeCpctDatamodel(final Iterable<EcrfDatamodelField> datamodelFields)
    {
        writeDatamodel(datamodelFields, this::cpctDatamodelInserter);
    }

    void writeDrupDatamodel(final Iterable<EcrfDatamodelField> datamodelFields)
    {
        writeDatamodel(datamodelFields, this::drupDatamodelInserter);
    }

    private void writeDatamodel(final Iterable<EcrfDatamodelField> datamodelFields, final Supplier<InsertValuesStep4> inserter)
    {
        context.batch(StreamSupport.stream(datamodelFields.spliterator(), false).map(field ->
        {
            String codeList = StringUtils.join(field.codeList().values().toArray(), ",");
            return inserter.get().values(field.name(), field.description(), codeList, field.isRelevant() ? "TRUE" : "FALSE");
        }).collect(Collectors.toList())).execute();
    }

    void writeCpctPatient(final EcrfPatient patient, boolean sequenced)
    {
        writePatient(patient, sequenced, this::cpctInserter);
    }

    void writeDrupPatient(final EcrfPatient patient, boolean sequenced)
    {
        writePatient(patient, sequenced, this::drupInserter);
    }

    private void writePatient(final EcrfPatient patient, boolean sequenced, final Supplier<InsertValuesStep14> inserter)
    {
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

    private InsertValuesStep4 cpctDatamodelInserter()
    {
        return context.insertInto(CPCTECRFDATAMODEL,
                CPCTECRFDATAMODEL.FIELDNAME,
                CPCTECRFDATAMODEL.DESCRIPTION,
                CPCTECRFDATAMODEL.CODELIST,
                CPCTECRFDATAMODEL.RELEVANT);
    }

    private InsertValuesStep4 drupDatamodelInserter()
    {
        return context.insertInto(DRUPECRFDATAMODEL,
                DRUPECRFDATAMODEL.FIELDNAME,
                DRUPECRFDATAMODEL.DESCRIPTION,
                DRUPECRFDATAMODEL.CODELIST,
                DRUPECRFDATAMODEL.RELEVANT);
    }

    private InsertValuesStep14 cpctInserter()
    {
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

    private InsertValuesStep14 drupInserter()
    {
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
