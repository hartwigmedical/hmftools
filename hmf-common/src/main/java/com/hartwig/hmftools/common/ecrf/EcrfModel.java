package com.hartwig.hmftools.common.ecrf;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.List;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfDatamodelField;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusModel;
import com.hartwig.hmftools.common.ecrf.formstatus.ImmutableFormStatusModel;
import com.hartwig.hmftools.common.ecrf.reader.XMLEcrfDatamodel;
import com.hartwig.hmftools.common.ecrf.reader.XMLEcrfDatamodelReader;
import com.hartwig.hmftools.common.ecrf.reader.XMLEcrfDatamodelToEcrfFields;
import com.hartwig.hmftools.common.ecrf.reader.XMLPatientReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EcrfModel {

    private static final Logger LOGGER = LogManager.getLogger(EcrfModel.class);

    @NotNull
    private final XMLEcrfDatamodel datamodel;
    @NotNull
    private final Iterable<EcrfDatamodelField> fields;
    @NotNull
    private final Iterable<EcrfPatient> patients;

    @NotNull
    public static EcrfModel loadFromXML(@NotNull final String ecrfXmlPath) throws XMLStreamException, FileNotFoundException {
        return loadFromXMLWithFormStates(ecrfXmlPath, new ImmutableFormStatusModel(Maps.newHashMap()));
    }

    @NotNull
    public static EcrfModel loadFromXMLWithFormStates(@NotNull final String ecrfXmlPath, @NotNull final FormStatusModel formStatusModel)
            throws XMLStreamException, FileNotFoundException {
        final XMLInputFactory factory = XMLInputFactory.newInstance();
        final XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(ecrfXmlPath));
        final XMLEcrfDatamodel datamodel = XMLEcrfDatamodelReader.readXMLDatamodel(reader);
        final Iterable<EcrfPatient> patients = XMLPatientReader.readPatients(reader, datamodel, formStatusModel);

        return new EcrfModel(datamodel, patients);
    }

    private EcrfModel(@NotNull final XMLEcrfDatamodel datamodel, @NotNull final Iterable<EcrfPatient> patients) {
        this.datamodel = datamodel;
        this.patients = patients;
        this.fields = XMLEcrfDatamodelToEcrfFields.convert(datamodel);
    }

    @NotNull
    public Iterable<EcrfDatamodelField> fields() {
        return fields;
    }

    public int patientCount() {
        return Lists.newArrayList(patients).size();
    }

    @NotNull
    public Iterable<EcrfPatient> findPatientsById(@NotNull final Iterable<String> patientIds) {
        final List<EcrfPatient> filteredPatients = Lists.newArrayList();
        for (final String patientId : patientIds) {
            final EcrfPatient patient = findPatientById(patientId);
            if (patient != null) {
                filteredPatients.add(patient);
            } else {
                LOGGER.warn("Did not find patient " + patientId + ": Adding dummy patient.");
                filteredPatients.add(new EcrfPatient(patientId, Maps.newHashMap(), Lists.newArrayList()));
            }
        }
        return filteredPatients;
    }

    @Nullable
    private EcrfPatient findPatientById(@NotNull final String patientId) {
        for (final EcrfPatient patient : patients) {
            if (patient.patientId().equals(patientId)) {
                return patient;
            }
        }
        return null;
    }

    @NotNull
    public Iterable<EcrfDatamodelField> findFieldsById(@NotNull final Iterable<String> fieldIds) {
        final List<EcrfDatamodelField> filteredFields = Lists.newArrayList();
        for (final String fieldId : fieldIds) {
            final EcrfDatamodelField field = findFieldById(fieldId);
            if (field != null) {
                filteredFields.add(field);
            } else {
                LOGGER.warn("Did not find field " + fieldId);
            }
        }
        return filteredFields;
    }

    @Nullable
    private EcrfDatamodelField findFieldById(@NotNull final String fieldId) {
        for (final EcrfDatamodelField field : fields) {
            if (field.name().equals(fieldId)) {
                return field;
            }
        }
        return null;
    }

    @NotNull
    public Iterable<EcrfPatient> patients() {
        return patients;
    }

    @NotNull
    public XMLEcrfDatamodel datamodel() {
        return datamodel;
    }
}
