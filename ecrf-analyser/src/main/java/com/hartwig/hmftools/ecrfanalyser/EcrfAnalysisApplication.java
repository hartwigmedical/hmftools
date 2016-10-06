package com.hartwig.hmftools.ecrfanalyser;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfField;
import com.hartwig.hmftools.ecrfanalyser.reader.EcrfFieldReader;

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
        //        EcrfPatientReader.readPatients(reader, allFields);

        BufferedWriter writer = new BufferedWriter(new FileWriter(csvOutPath, false));
        writer.write("CATEGORY, FIELDNAME, DESCRIPTION, VALUES");

        for (EcrfField field : allFields) {
            writer.newLine();
            writer.write(toCSV(field));
        }
        writer.close();
    }

    @NotNull
    private static String toCSV(@NotNull EcrfField field) {
        String valuesString = "";
        for (String value : field.values().values()) {
            valuesString += value.replaceAll(",", ":") + "; ";
        }
        if (valuesString.length() > 0) {
            valuesString = valuesString.substring(0, valuesString.length() - 2);
        }
        return field.category() + ", " + field.fieldName() + ", " + field.description().replaceAll(",", ":") + ", "
                + valuesString;
    }
}
