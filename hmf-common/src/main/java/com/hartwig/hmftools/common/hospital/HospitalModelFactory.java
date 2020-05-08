package com.hartwig.hmftools.common.hospital;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HospitalModelFactory {

    private static final Logger LOGGER = LogManager.getLogger(HospitalModelFactory.class);

    private static final String SAMPLE_HOSPITAL_MAPPING_CSV = "sample_hospital_mapping.csv";
    private static final String HOSPITALS_ADRESS_CSV = "hospital_adress.csv";
    private static final String HOSPITAL_CPCT_CSV = "hospital_cpct.csv";
    private static final String HOSPITAL_DRUP_CSV = "hospital_drup.csv";
    private static final String HOSPITAL_WIDE_CSV = "hospital_wide.csv";

    private static final int SAMPLE_MAPPING_ID_COLUMN = 0;
    private static final int HOSPITAL_MAPPING_COLUMN = 1;
    private static final int FIELD_COUNT_SAMPLE_HOSPITAL_MAPPING = 2;

    private static final int HOSPITAL_ADRESS_ID_COLUMN = 0;
    private static final int HOSPITAL_ADRESS_NAME_COLUMN = 1;
    private static final int HOSPITAL_ADRESS_ZIP_COLUMN = 2;
    private static final int HOSPITAL_ADRESS_CITY_COLUMN = 3;
    private static final int HOSPITAL_ADRESS_FIELD_COUNT = 4;

    private static final int HOSPITAL_DATA_ID_COLUMN = 0;
    private static final int HOSPITAL_DATA_PI_COLUMN = 1;
    private static final int HOSPITAL_DATA_REQUEST_NAME_COLUMN = 2;
    private static final int HOSPITAL_DATA_REQUEST_EMAIL_COLUMN = 3;
    private static final int HOSPITAL_DATA_FIELD_COUNT = 4;

    private static final String FIELD_SEPARATOR = ",";

    private HospitalModelFactory() {
    }

    @NotNull
    public static HospitalModel fromHospitalDirectory(@NotNull String hospitalDirectory) throws IOException {
        String sampleHospitalMappingCsv = hospitalDirectory + File.separator + SAMPLE_HOSPITAL_MAPPING_CSV;
        String hospitalDataCPCTCsv = hospitalDirectory + File.separator + HOSPITAL_CPCT_CSV;
        String hospitalDataDRUPCsv = hospitalDirectory + File.separator + HOSPITAL_DRUP_CSV;
        String hospitalDataWIDECsv = hospitalDirectory + File.separator + HOSPITAL_WIDE_CSV;

        String hospitalAdressCsv = hospitalDirectory + File.separator + HOSPITALS_ADRESS_CSV;

        Map<String, HospitalSampleMapping> sampleHospitalMapping = readFromSampleHospitalMapping(sampleHospitalMappingCsv);

        Map<String, HospitalData> hospitalDataCPCT = readFromHospitalDataCPCT(hospitalDataCPCTCsv);
        Map<String, HospitalData> hospitalDataDRUP = readFromHospitalDataDRUP(hospitalDataDRUPCsv);
        Map<String, HospitalData> hospitalDataWIDE = readFromHospitalDataWIDE(hospitalDataWIDECsv);

        Map<String, HospitalAdress> hospitalAdress = readFromHospitalAdress(hospitalAdressCsv);

        return ImmutableHospitalModel.of(
                sampleHospitalMapping,
                hospitalDataCPCT,
                hospitalDataDRUP,
                hospitalDataWIDE,
                hospitalAdress);
    }

    @NotNull
    public static HospitalModel empty() {
        return ImmutableHospitalModel.of(Maps.newHashMap(),
                Maps.newHashMap(),
                Maps.newHashMap(),
                Maps.newHashMap(),
                Maps.newHashMap());
    }


    @NotNull
    private static Map<String, HospitalSampleMapping> readFromSampleHospitalMapping(@NotNull String sampleHospitalMappingCsv)
            throws IOException {
        Map<String, HospitalSampleMapping> hospitalPerSampleMap = Maps.newHashMap();

        List<String> lines = FileReader.build().readLines(new File(sampleHospitalMappingCsv).toPath());
        for (String line : lines) {
            String[] parts = line.split(FIELD_SEPARATOR, FIELD_COUNT_SAMPLE_HOSPITAL_MAPPING);
            if (parts.length == FIELD_COUNT_SAMPLE_HOSPITAL_MAPPING) {
                HospitalSampleMapping hospitalManual = ImmutableHospitalSampleMapping.of(parts[HOSPITAL_MAPPING_COLUMN]);

                hospitalPerSampleMap.put(parts[SAMPLE_MAPPING_ID_COLUMN], hospitalManual);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in sample hospital mapping csv: '{}'", line);
            }
        }
        return hospitalPerSampleMap;
    }

    @NotNull
    private static Map<String, HospitalData> readFromHospitalDataCPCT(@NotNull String hospitalDataCsv) throws IOException {
        Map<String, HospitalData> hospitalData = Maps.newHashMap();
        List<String> lines = FileReader.build().readLines(new File(hospitalDataCsv).toPath());
        for (String line : lines) {
            String[] parts = line.split(FIELD_SEPARATOR, HOSPITAL_DATA_FIELD_COUNT);
            if (parts.length == HOSPITAL_DATA_FIELD_COUNT) {
                HospitalData hospital = ImmutableHospitalData.of(parts[HOSPITAL_DATA_ID_COLUMN],
                        parts[HOSPITAL_DATA_PI_COLUMN],
                        parts[HOSPITAL_DATA_REQUEST_NAME_COLUMN],
                        parts[HOSPITAL_DATA_REQUEST_EMAIL_COLUMN]);

                hospitalData.put(parts[HOSPITAL_DATA_ID_COLUMN], hospital);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in hospital csv: '{}'", line);
            }
        }
        return hospitalData;
    }

    @NotNull
    private static Map<String, HospitalData> readFromHospitalDataDRUP(@NotNull String hospitalDataCsv) throws IOException {
        Map<String, HospitalData> hospitalData = Maps.newHashMap();
        List<String> lines = FileReader.build().readLines(new File(hospitalDataCsv).toPath());
        for (String line : lines) {
            String[] parts = line.split(FIELD_SEPARATOR, HOSPITAL_DATA_FIELD_COUNT);
            if (parts.length == HOSPITAL_DATA_FIELD_COUNT) {
                HospitalData hospital = ImmutableHospitalData.of(parts[HOSPITAL_DATA_ID_COLUMN],
                        parts[HOSPITAL_DATA_PI_COLUMN],
                        parts[HOSPITAL_DATA_REQUEST_NAME_COLUMN],
                        parts[HOSPITAL_DATA_REQUEST_EMAIL_COLUMN]);

                hospitalData.put(parts[HOSPITAL_DATA_ID_COLUMN], hospital);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in hospital csv: '{}'", line);
            }
        }
        return hospitalData;
    }

    @NotNull
    private static Map<String, HospitalData> readFromHospitalDataWIDE(@NotNull String hospitalDataCsv) throws IOException {
        Map<String, HospitalData> hospitalData = Maps.newHashMap();
        List<String> lines = FileReader.build().readLines(new File(hospitalDataCsv).toPath());
        for (String line : lines) {
            String[] parts = line.split(FIELD_SEPARATOR, HOSPITAL_DATA_FIELD_COUNT);
            if (parts.length == HOSPITAL_DATA_FIELD_COUNT) {
                HospitalData hospital = ImmutableHospitalData.of(parts[HOSPITAL_DATA_ID_COLUMN],
                        parts[HOSPITAL_DATA_PI_COLUMN],
                        parts[HOSPITAL_DATA_REQUEST_NAME_COLUMN],
                        parts[HOSPITAL_DATA_REQUEST_EMAIL_COLUMN]);

                hospitalData.put(parts[HOSPITAL_DATA_ID_COLUMN], hospital);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in hospital csv: '{}'", line);
            }
        }
        return hospitalData;
    }

    @NotNull
    private static Map<String, HospitalAdress> readFromHospitalAdress(@NotNull String hospitalAdressCsv) throws IOException {
        Map<String, HospitalAdress> hospitalAdress = Maps.newHashMap();
        List<String> lines = FileReader.build().readLines(new File(hospitalAdressCsv).toPath());
        for (String line : lines) {
            String[] parts = line.split(FIELD_SEPARATOR, HOSPITAL_ADRESS_FIELD_COUNT);
            if (parts.length == HOSPITAL_ADRESS_FIELD_COUNT) {
                HospitalAdress hospital = ImmutableHospitalAdress.of(parts[HOSPITAL_ADRESS_ID_COLUMN],
                        parts[HOSPITAL_ADRESS_NAME_COLUMN],
                        parts[HOSPITAL_ADRESS_ZIP_COLUMN],
                        parts[HOSPITAL_ADRESS_CITY_COLUMN]);

                hospitalAdress.put(parts[SAMPLE_MAPPING_ID_COLUMN], hospital);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in hospital adress csv: '{}'", line);
            }
        }
        return hospitalAdress;
    }
}
