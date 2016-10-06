package com.hartwig.hmftools.ecrfanalyser;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfField;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfPatient;
import com.hartwig.hmftools.ecrfanalyser.reader.EcrfFieldReader;
import com.hartwig.hmftools.ecrfanalyser.reader.EcrfPatientReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class EcrfAnalysisApplication {

    private static final Logger LOGGER = LogManager.getLogger(EcrfAnalysisApplication.class);

    public static void main(final String... args) {
        LOGGER.info("Hello World");
    }

    @NotNull
    private final String ecrfXmlPath;
    @NotNull
    private final String csvOutPath;

    EcrfAnalysisApplication(@NotNull final String ecrfXmlPath, @NotNull final String csvOutPath) {
        this.ecrfXmlPath = ecrfXmlPath;
        this.csvOutPath = csvOutPath;
    }

    void generateCsv(@NotNull List<String> patients)
            throws IOException, XMLStreamException {
        XMLInputFactory factory = XMLInputFactory.newInstance();
        XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(ecrfXmlPath));

        List<EcrfField> allFields = EcrfFieldReader.readFields(reader);
        List<EcrfPatient> allPatients = EcrfPatientReader.readPatients(reader, allFields);

        writePatientsToCSV(allPatients, csvOutPath);
    }

    private static void writePatientsToCSV(@NotNull List<EcrfPatient> patients, @NotNull String csvOutPath)
            throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(csvOutPath, false));

        Set<EcrfField> fields = Sets.newTreeSet(patients.get(0).fields());

        String header = "PATIENT";

        for (EcrfPatient patient : patients) {
            header += (", " + patient.patientId());
        }
        writer.write(header);

        for (EcrfField field : fields) {
            writer.newLine();
            // KODU: Pass fields of first patient to enforce consistency across all patients.
            writer.write(fieldDataToCSV(field, patients));
        }
        writer.close();
    }

    @NotNull
    private static String fieldDataToCSV(@NotNull EcrfField field, @NotNull List<EcrfPatient> patients) {
        String fieldCSV = toFieldName(field);
        for (EcrfPatient patient : patients) {
            String value = patient.fieldValue(field);
            fieldCSV += ", " + value.replaceAll(",", ":");
        }
        return fieldCSV;
    }

    @NotNull
    private static String toFieldName(@NotNull EcrfField field) {
        return field.category() + "." + field.fieldName();
    }

    @SuppressWarnings("unused")
    private static void writeDatamodelToCSV(@NotNull List<EcrfField> fields, @NotNull String csvOutPath)
            throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(csvOutPath, false));
        writer.write("CATEGORY, FIELDNAME, DESCRIPTION, VALUES");

        for (EcrfField field : fields) {
            writer.newLine();
            writer.write(fieldToCSV(field));
        }
        writer.close();
    }

    @NotNull
    private static String fieldToCSV(@NotNull EcrfField field) {
        String valuesString = "";
        for (String value : field.values().values()) {
            valuesString += value.replaceAll(",", ":") + ", ";
        }
        if (valuesString.length() > 0) {
            valuesString = valuesString.substring(0, valuesString.length() - 2);
        }
        return field.category() + ", " + field.fieldName() + ", " + field.description().replaceAll(",", ":") + ", "
                + valuesString;
    }
}
