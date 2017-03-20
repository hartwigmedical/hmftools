package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfField;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class PatientDbRunner {
    private static final Logger LOGGER = LogManager.getLogger(PatientDbRunner.class);
    private static final String RUNS_DIR = "runsDir";
    private static final String ECRF_FILE = "ecrf";

    public static void main(String[] args)
            throws ParseException, IOException, InterruptedException, java.text.ParseException, XMLStreamException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String runFolderPath = cmd.getOptionValue(RUNS_DIR);
        final String ecrfFilePath = cmd.getOptionValue(ECRF_FILE);

        if (runFolderPath == null) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-Db", options);
        } else {
            final File dir = new File(runFolderPath);
            if (dir.isDirectory()) {
                final List<CpctRunData> data = RunsFolderProcessor.getPatientRunsData(dir);
                LOGGER.info("Listing data for " + data.size() + " patients.");
                LOGGER.info(data.toString());
                final DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");
                LOGGER.info("Loading ecrf model...");
                final CpctEcrfModel model = CpctEcrfModel.loadFromXML(ecrfFilePath);
                final List<String> cpctPatientIds = data.stream().map(CpctRunData::getPatientId).filter(
                        id -> id.startsWith("CPCT")).collect(Collectors.toList());
                final Iterable<EcrfPatient> patients = model.findPatientsById(cpctPatientIds);
                final Map<Integer, String> hospitals = getHospitals(model);
                LOGGER.info("Reading patient data...");
                for (final EcrfPatient patient : patients) {
                    final String sex = getField(patient, "BASELINE.DEMOGRAPHY.DEMOGRAPHY.SEX");
                    final String ethnicity = getField(patient, "BASELINE.DEMOGRAPHY.DEMOGRAPHY.ETHNIC");
                    final Integer birthYear = getBirthYear(patient, dateFormat);
                    final String hospital = getHospital(patient, hospitals);
                    final PatientData patientData = new PatientData(patient.patientId(), null, sex, birthYear,
                            hospital, ethnicity);
                    final String tumorLocation = getField(patient, "BASELINE.CARCINOMA.CARCINOMA.PTUMLOC");
                    final String tumorEntryStage = getField(patient, "BASELINE.CARCINOMA.CARCINOMA.ENTRYSTAGE");
                    final List<String> biopsyLocations = getFieldValues(patient, "BIOPSY.BIOPS.BIOPSIES.BILESSITE");
                    final TumorData tumorData = new TumorData(tumorLocation, biopsyLocations, tumorEntryStage);
                    LOGGER.info(patientData.toString());
                    LOGGER.info(tumorData.toString());
                }
            } else {
                if (!dir.exists()) {
                    LOGGER.warn("dir " + dir + " does not exist.");
                }
                HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp("Patient-Db", options);
            }
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(RUNS_DIR, true, "Path towards the folder containing patient runs.");
        options.addOption(ECRF_FILE, true, "Path towards the cpct ecrf file.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @Nullable
    private static String getField(@NotNull EcrfPatient patient, @NotNull String fieldName) {
        final List<String> values = patient.fieldValuesByName(fieldName);
        if (values == null) {
            LOGGER.warn(fieldName + " not found for patient " + patient.patientId() + ".");
            return null;
        } else if (values.size() == 0) {
            LOGGER.warn(fieldName + " for patient " + patient.patientId() + " contains no values.");
            return null;
        } else if (values.size() > 1) {
            LOGGER.warn(fieldName + " for patient " + patient.patientId() + " contains more than 1 value.");
        }
        return values.get(0);
    }

    @NotNull
    private static List<String> getFieldValues(@NotNull EcrfPatient patient, @NotNull String fieldName) {
        final List<String> values = patient.fieldValuesByName(fieldName);
        if (values == null) {
            LOGGER.warn(fieldName + " not found for patient " + patient.patientId() + ".");
            return null;
        } else if (values.size() == 0) {
            LOGGER.warn(fieldName + " for patient " + patient.patientId() + " contains no values.");
            return null;
        } else
            return values;
    }

    @Nullable
    private static Integer getBirthYear(@NotNull EcrfPatient patient, @NotNull DateFormat dateFormat) {
        final String nbirthYear = getField(patient, "BASELINE.SELCRIT.SELCRIT.NBIRTHYEAR");
        if (nbirthYear != null && !nbirthYear.equals("")) {
            return Integer.parseInt(nbirthYear);
        }
        final String birthYear = getField(patient, "BASELINE.ELIGIBILITY.ELIGIBILITY.BIRTHYEAR");
        if (birthYear != null && !birthYear.equals("")) {
            return Integer.parseInt(birthYear);
        }
        final String birthdtces = getField(patient, "BASELINE.ELIGIBILITY.ELIGIBILITY.BIRTHDTCES");
        if (birthdtces != null && !birthdtces.equals("")) {
            try {
                final Date date = dateFormat.parse(birthdtces);
                final Calendar calendar = Calendar.getInstance();
                calendar.setTime(date);
                return calendar.get(Calendar.YEAR);
            } catch (java.text.ParseException e) {
                LOGGER.info("BASELINE.ELIGIBILITY.ELIGIBILITY.BIRTHDTCES field did not contain valid date for patient "
                        + patient.patientId());
                return null;
            }
        }
        LOGGER.info("No completed birth year fields found for patient " + patient.patientId());
        return null;
    }

    @NotNull
    private static String getHospital(@NotNull EcrfPatient patient, @NotNull Map<Integer, String> hospitals) {
        final Integer hospitalCode = Integer.parseInt(patient.patientId().substring(6, 8));
        return hospitals.get(hospitalCode);
    }

    @NotNull
    private static Map<Integer, String> getHospitals(@NotNull CpctEcrfModel datamodel) {
        final Map<Integer, String> hospitals = Maps.newHashMap();
        final Iterable<EcrfField> fields = datamodel.findFieldsById(
                Lists.newArrayList("BASELINE.ELIGIBILITY.ELIGIBILITY.HOSPITAL", "BASELINE.SELCRIT.SELCRIT.NHOSPITAL"));
        fields.forEach(field -> hospitals.putAll(field.codeList()));
        return hospitals;
    }
}