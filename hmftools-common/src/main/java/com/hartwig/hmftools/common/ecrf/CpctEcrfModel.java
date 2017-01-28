package com.hartwig.hmftools.common.ecrf;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.List;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfField;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.reader.XMLEcrfDatamodel;
import com.hartwig.hmftools.common.ecrf.reader.XMLEcrfDatamodelReader;
import com.hartwig.hmftools.common.ecrf.reader.XMLEcrfDatamodelToEcrfFields;
import com.hartwig.hmftools.common.ecrf.reader.XMLPatientReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CpctEcrfModel {

    private static final Logger LOGGER = LogManager.getLogger(CpctEcrfModel.class);

    @NotNull
    private final Iterable<EcrfField> fields;
    @NotNull
    private final Iterable<EcrfPatient> patients;

    @NotNull
    public static CpctEcrfModel loadFromXML(@NotNull final String ecrfXmlPath)
            throws XMLStreamException, FileNotFoundException {
        final XMLInputFactory factory = XMLInputFactory.newInstance();
        final XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(ecrfXmlPath));
        final XMLEcrfDatamodel datamodel = XMLEcrfDatamodelReader.readXMLDatamodel(reader);

        final Iterable<EcrfField> fields = Sets.newTreeSet(XMLEcrfDatamodelToEcrfFields.convert(datamodel));
        final Iterable<EcrfPatient> patients = XMLPatientReader.readPatients(reader, fields);

        return new CpctEcrfModel(fields, patients);
    }

    public CpctEcrfModel(@NotNull final Iterable<EcrfField> fields, @NotNull final Iterable<EcrfPatient> patients) {
        this.fields = fields;
        this.patients = patients;
    }

    @NotNull
    public Iterable<EcrfField> fields() {
        return fields;
    }

    @NotNull
    public Iterable<EcrfPatient> findPatientsById(@NotNull Iterable<String> patientIds) {
        final List<EcrfPatient> filteredPatients = Lists.newArrayList();
        for (final String patientId : patientIds) {
            final EcrfPatient patient = findPatientById(patientId);
            if (patient != null) {
                filteredPatients.add(patient);
            } else {
                LOGGER.warn("Did not find patient " + patientId + ": Adding dummy patient.");
                filteredPatients.add(new EcrfPatient(patientId, Maps.newHashMap()));
            }
        }
        return filteredPatients;

    }
    @Nullable
    public EcrfPatient findPatientById(@NotNull final String patientId) {
        for (final EcrfPatient patient : patients) {
            if (patient.patientId().equals(patientId)) {
                return patient;
            }
        }
        return null;
    }

    @NotNull
    public Iterable<EcrfField> findFieldsById(@NotNull Iterable<String> fieldIds) {
        final List<EcrfField> filteredFields = Lists.newArrayList();
        for (final String fieldId : fieldIds) {
            final EcrfField field = findFieldById(fieldId);
            if (field != null) {
                filteredFields.add(field);
            } else {
                LOGGER.warn("Did not find field " + fieldId);
            }
        }
        return filteredFields;
    }

    @Nullable
    public EcrfField findFieldById(@NotNull final String fieldId) {
        for (final EcrfField field : fields) {
            if (field.name().equals(fieldId)) {
                return field;
            }
        }
        return null;
    }
}
