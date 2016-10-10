package com.hartwig.hmftools.ecrfanalyser;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfField;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfPatient;
import com.hartwig.hmftools.ecrfanalyser.reader.XMLEcrfDatamodel;
import com.hartwig.hmftools.ecrfanalyser.reader.XMLEcrfDatamodelReader;
import com.hartwig.hmftools.ecrfanalyser.reader.XMLEcrfDatamodelToEcrfFields;
import com.hartwig.hmftools.ecrfanalyser.reader.XMLPatientReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EcrfAnalysisApplication {

    private static final Logger LOGGER = LogManager.getLogger(EcrfAnalysisApplication.class);

    private static final String ECRF_XML_PATH_ARGS_DESC = "The path to the ecrf xml file";
    private static final String ECRF_XML_PATH = "ecrf";

    private static final String CSV_OUT_PATH_ARGS_DESC = "The file to write csv results to.";
    private static final String CSV_OUT_PATH = "csvout";

    private static final String PATIENTS_ARGS_DESC = "A comma-separated list of patients to filter on";
    private static final String PATIENTS = "patients";

    public static void main(final String... args) throws ParseException, IOException, XMLStreamException {
        Options options = createOptions();
        CommandLine cmd = createCommandLine(options, args);

        String reportsPath = cmd.getOptionValue(ECRF_XML_PATH);
        String csvOut = cmd.getOptionValue(CSV_OUT_PATH);
        String patients = cmd.getOptionValue(PATIENTS);

        if (reportsPath == null || csvOut == null || patients == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Ecrf-Analyser", options);
            System.exit(1);
        }

        new EcrfAnalysisApplication(reportsPath, csvOut, Lists.newArrayList(patients.split(","))).run();
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(ECRF_XML_PATH, true, ECRF_XML_PATH_ARGS_DESC);
        options.addOption(CSV_OUT_PATH, true, CSV_OUT_PATH_ARGS_DESC);
        options.addOption(PATIENTS, true, PATIENTS_ARGS_DESC);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private final String ecrfXmlPath;
    @NotNull
    private final String csvOutPath;
    @NotNull
    private final List<String> patientIds;

    EcrfAnalysisApplication(@NotNull final String ecrfXmlPath, @NotNull final String csvOutPath,
            @NotNull final List<String> patientIds) {
        this.ecrfXmlPath = ecrfXmlPath;
        this.csvOutPath = csvOutPath;
        this.patientIds = patientIds;
    }

    void run() throws IOException, XMLStreamException {
        XMLInputFactory factory = XMLInputFactory.newInstance();
        XMLStreamReader reader = factory.createXMLStreamReader(new FileInputStream(ecrfXmlPath));
        XMLEcrfDatamodel datamodel = XMLEcrfDatamodelReader.readXMLDatamodel(reader);

        Iterable<EcrfField> allFields = Sets.newTreeSet(XMLEcrfDatamodelToEcrfFields.convert(datamodel));
        Iterable<EcrfPatient> allPatients = XMLPatientReader.readPatients(reader, allFields);

        List<EcrfPatient> filteredPatients = Lists.newArrayList();
        for (String patientId : patientIds) {
            EcrfPatient patient = findPatient(allPatients, patientId);
            if (patient != null) {
                filteredPatients.add(patient);
            } else {
                LOGGER.warn("Did not find patient " + patientId);
                filteredPatients.add(new EcrfPatient(patientId, Maps.<EcrfField, List<String>>newHashMap()));
            }
        }

        //        writeDatamodelToCSV(allFields, csvOutPath);
        writePatientsToCSV(filteredPatients, allFields, csvOutPath);
    }

    @Nullable
    private EcrfPatient findPatient(@NotNull Iterable<EcrfPatient> patients, @NotNull String patientIdToFind) {
        for (EcrfPatient patient : patients) {
            if (patient.patientId().equals(patientIdToFind)) {
                return patient;
            }
        }
        return null;
    }

    @SuppressWarnings("unused")
    private static void writePatientsToCSV(@NotNull List<EcrfPatient> patients, @NotNull Iterable<EcrfField> fields,
            @NotNull String csvOutPath) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(csvOutPath, false));

        String header = "PATIENT";

        for (EcrfPatient patient : patients) {
            header += (", " + patient.patientId());
        }
        writer.write(header);

        for (EcrfField field : fields) {
            writer.newLine();
            writer.write(fieldDataToCSV(field, patients));
        }
        writer.close();
    }

    @NotNull
    private static String fieldDataToCSV(@NotNull EcrfField field, @NotNull List<EcrfPatient> patients) {
        String fieldCSV = field.name();
        for (EcrfPatient patient : patients) {
            List<String> values = patient.fieldValues(field);
            String finalValue = Strings.EMPTY;
            if (values != null && containsSomeValue(values)) {
                for (int i = 0; i < values.size(); i++) {
                    if (i == 0) {
                        finalValue = values.get(i);
                    } else {
                        finalValue += ("; " + values.get(i));
                    }
                }
            }
            fieldCSV += ", " + finalValue.replaceAll(",", ":");
        }
        return fieldCSV;
    }

    private static boolean containsSomeValue(@NotNull Iterable<String> strings) {
        for (String value : strings) {
            if (!value.isEmpty()) {
                return true;
            }
        }
        return false;
    }

    @SuppressWarnings("unused")
    private static void writeDatamodelToCSV(@NotNull Iterable<EcrfField> fields, @NotNull String csvOutPath)
            throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(csvOutPath, false));
        writer.write("FIELD, DESCRIPTION, VALUES");

        for (EcrfField field : fields) {
            writer.newLine();
            writer.write(fieldToCSV(field));
        }
        writer.close();
    }

    @NotNull
    private static String fieldToCSV(@NotNull EcrfField field) {
        String valuesString = "";
        for (String value : field.codeList().values()) {
            valuesString += value.replaceAll(",", ":") + "; ";
        }
        if (valuesString.length() > 0) {
            valuesString = valuesString.substring(0, valuesString.length() - 2);
        }
        return field.name() + ", " + field.description().replaceAll(",", ":") + ", " + valuesString;
    }
}
