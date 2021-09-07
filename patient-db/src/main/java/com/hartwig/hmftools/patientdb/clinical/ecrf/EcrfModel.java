package com.hartwig.hmftools.patientdb.clinical.ecrf;

import java.io.FileInputStream;
import java.io.FileNotFoundException;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfDatamodelField;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatusModel;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.ImmutableFormStatusModel;
import com.hartwig.hmftools.patientdb.clinical.ecrf.reader.XMLEcrfDatamodel;
import com.hartwig.hmftools.patientdb.clinical.ecrf.reader.XMLEcrfDatamodelReader;
import com.hartwig.hmftools.patientdb.clinical.ecrf.reader.XMLEcrfDatamodelToEcrfFields;
import com.hartwig.hmftools.patientdb.clinical.ecrf.reader.XMLPatientReader;

import org.jetbrains.annotations.NotNull;

public class EcrfModel {

    @NotNull
    private final XMLEcrfDatamodel datamodel;
    @NotNull
    private final Iterable<EcrfDatamodelField> fields;
    @NotNull
    private final Iterable<EcrfPatient> patients;

    @NotNull
    public static EcrfModel loadFromXMLNoFormStates(@NotNull String ecrfXmlPath) throws XMLStreamException, FileNotFoundException {
        return loadFromXMLWithFormStates(ecrfXmlPath, new ImmutableFormStatusModel(Maps.newHashMap()));
    }

    @NotNull
    public static EcrfModel loadFromXMLWithFormStates(@NotNull String ecrfXmlPath, @NotNull FormStatusModel formStatusModel)
            throws XMLStreamException, FileNotFoundException {
        XMLInputFactory factory = XMLInputFactory.newInstance();
        XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(ecrfXmlPath));
        XMLEcrfDatamodel datamodel = XMLEcrfDatamodelReader.readXMLDatamodel(reader);
        Iterable<EcrfPatient> patients = XMLPatientReader.readPatients(reader, datamodel, formStatusModel);

        return new EcrfModel(datamodel, patients);
    }

    private EcrfModel(@NotNull final XMLEcrfDatamodel datamodel, @NotNull final Iterable<EcrfPatient> patients) {
        this.datamodel = datamodel;
        this.patients = patients;
        this.fields = XMLEcrfDatamodelToEcrfFields.convert(datamodel);
    }

    @NotNull
    public Iterable<EcrfPatient> patients() {
        return patients;
    }

    public int patientCount() {
        return Lists.newArrayList(patients).size();
    }

    @NotNull
    public Iterable<EcrfDatamodelField> fields() {
        return fields;
    }

    @NotNull
    public XMLEcrfDatamodel datamodel() {
        return datamodel;
    }
}
