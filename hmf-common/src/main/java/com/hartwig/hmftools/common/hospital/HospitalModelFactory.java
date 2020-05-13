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

    private static final String SAMPLE_HOSPITAL_MAPPING_CSV = "sample_hospital_mapping.tsv";
    private static final String HOSPITALS_ADDRESS_CSV = "hospital_address.tsv";
    private static final String HOSPITAL_CPCT_CSV = "hospital_cpct.tsv";
    private static final String HOSPITAL_DRUP_CSV = "hospital_drup.tsv";
    private static final String HOSPITAL_WIDE_CSV = "hospital_wide.tsv";

    private static final int SAMPLE_MAPPING_ID_COLUMN = 0;
    private static final int HOSPITAL_MAPPING_COLUMN = 1;
    private static final int FIELD_COUNT_SAMPLE_HOSPITAL_MAPPING = 2;

    private static final int HOSPITAL_ADDRESS_ID_COLUMN = 0;
    private static final int HOSPITAL_ADDRESS_NAME_COLUMN = 1;
    private static final int HOSPITAL_ADDRESS_ZIP_COLUMN = 2;
    private static final int HOSPITAL_ADDRESS_CITY_COLUMN = 3;
    private static final int HOSPITAL_ADDRESS_FIELD_COUNT = 4;

    private static final int HOSPITAL_DATA_ID_COLUMN = 0;
    private static final int HOSPITAL_DATA_PI_COLUMN = 1;
    private static final int HOSPITAL_DATA_REQUEST_NAME_COLUMN = 2;
    private static final int HOSPITAL_DATA_REQUEST_EMAIL_COLUMN = 3;
    private static final int HOSPITAL_DATA_FIELD_COUNT_WIDE = 4;
    private static final int HOSPITAL_DATA_FIELD_COUNT_CPCT_DRUP = 2;

    private static final String FIELD_SEPARATOR = "\t";

    private HospitalModelFactory() {
    }

    @NotNull
    public static HospitalModel fromHospitalDirectory(@NotNull String limsDirectory) throws IOException {
        String sampleHospitalMappingCsv = limsDirectory + File.separator + SAMPLE_HOSPITAL_MAPPING_CSV;
        String hospitalDataCPCTCsv = limsDirectory + File.separator + HOSPITAL_CPCT_CSV;
        String hospitalDataDRUPCsv = limsDirectory + File.separator + HOSPITAL_DRUP_CSV;
        String hospitalDataWIDECsv = limsDirectory + File.separator + HOSPITAL_WIDE_CSV;

        String hospitalAddressCsv = limsDirectory + File.separator + HOSPITALS_ADDRESS_CSV;

        Map<String, HospitalSampleMapping> sampleHospitalMapping = readFromSampleHospitalMapping(sampleHospitalMappingCsv);

        Map<String, HospitalData> hospitalDataCPCT = readFromHospitalDataCPCT(hospitalDataCPCTCsv);
        Map<String, HospitalData> hospitalDataDRUP = readFromHospitalDataDRUP(hospitalDataDRUPCsv);
        Map<String, HospitalData> hospitalDataWIDE = readFromHospitalDataWIDE(hospitalDataWIDECsv);

        Map<String, HospitalAddress> hospitalAddress = readFromHospitalAddress(hospitalAddressCsv);

        return ImmutableHospitalModel.of(sampleHospitalMapping, hospitalDataCPCT, hospitalDataDRUP, hospitalDataWIDE, hospitalAddress);
    }

    @NotNull
    public static HospitalModel empty() {
        return ImmutableHospitalModel.of(Maps.newHashMap(), Maps.newHashMap(), Maps.newHashMap(), Maps.newHashMap(), Maps.newHashMap());
    }

    @NotNull
    public static Map<String, HospitalSampleMapping> readFromSampleHospitalMapping(@NotNull String sampleHospitalMappingCsv)
            throws IOException {
        Map<String, HospitalSampleMapping> hospitalPerSampleMap = Maps.newHashMap();

        List<String> lines = FileReader.build().readLines(new File(sampleHospitalMappingCsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_SEPARATOR, FIELD_COUNT_SAMPLE_HOSPITAL_MAPPING);
            if (parts.length == FIELD_COUNT_SAMPLE_HOSPITAL_MAPPING) {
                HospitalSampleMapping hospitalManual = ImmutableHospitalSampleMapping.of(parts[HOSPITAL_MAPPING_COLUMN]);

                hospitalPerSampleMap.put(parts[SAMPLE_MAPPING_ID_COLUMN], hospitalManual);
            } else {
                LOGGER.warn("Could not properly parse line in sample hospital mapping csv: '{}'", line);
            }
        }
        return hospitalPerSampleMap;
    }

    @NotNull
    public static Map<String, HospitalData> readFromHospitalDataCPCT(@NotNull String hospitalDataCsv) throws IOException {
        Map<String, HospitalData> hospitalData = Maps.newHashMap();
        List<String> lines = FileReader.build().readLines(new File(hospitalDataCsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_SEPARATOR, HOSPITAL_DATA_FIELD_COUNT_CPCT_DRUP);
            if (parts.length == HOSPITAL_DATA_FIELD_COUNT_CPCT_DRUP) {
                HospitalData hospital = ImmutableHospitalData.of(parts[HOSPITAL_DATA_ID_COLUMN],
                        parts[HOSPITAL_DATA_PI_COLUMN], null, null);

                hospitalData.put(parts[HOSPITAL_DATA_ID_COLUMN], hospital);
            } else {
                LOGGER.warn("Could not properly parse line in cpct hospital csv: '{}'", line);
            }
        }
        return hospitalData;
    }

    @NotNull
    public static Map<String, HospitalData> readFromHospitalDataDRUP(@NotNull String hospitalDataCsv) throws IOException {
        Map<String, HospitalData> hospitalData = Maps.newHashMap();
        List<String> lines = FileReader.build().readLines(new File(hospitalDataCsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_SEPARATOR, HOSPITAL_DATA_FIELD_COUNT_CPCT_DRUP);
            if (parts.length == HOSPITAL_DATA_FIELD_COUNT_CPCT_DRUP) {
                HospitalData hospital = ImmutableHospitalData.of(parts[HOSPITAL_DATA_ID_COLUMN],
                        parts[HOSPITAL_DATA_PI_COLUMN], null, null);

                hospitalData.put(parts[HOSPITAL_DATA_ID_COLUMN], hospital);
            } else {
                LOGGER.warn("Could not properly parse line in drup hospital csv: '{}'", line);
            }
        }
        return hospitalData;
    }

    @NotNull
    public static Map<String, HospitalData> readFromHospitalDataWIDE(@NotNull String hospitalDataCsv) throws IOException {
        Map<String, HospitalData> hospitalData = Maps.newHashMap();
        List<String> lines = FileReader.build().readLines(new File(hospitalDataCsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_SEPARATOR, HOSPITAL_DATA_FIELD_COUNT_WIDE);
            if (parts.length == HOSPITAL_DATA_FIELD_COUNT_WIDE) {
                HospitalData hospital = ImmutableHospitalData.of(parts[HOSPITAL_DATA_ID_COLUMN],
                        parts[HOSPITAL_DATA_PI_COLUMN],
                        parts[HOSPITAL_DATA_REQUEST_NAME_COLUMN],
                        parts[HOSPITAL_DATA_REQUEST_EMAIL_COLUMN]);

                hospitalData.put(parts[HOSPITAL_DATA_ID_COLUMN], hospital);
            } else {
                LOGGER.warn("Could not properly parse line in wide hospital csv: '{}'", line);
            }
        }
        return hospitalData;
    }

    @NotNull
    public static Map<String, HospitalAddress> readFromHospitalAddress(@NotNull String hospitalAddressCsv) throws IOException {
        Map<String, HospitalAddress> hospitalAddress = Maps.newHashMap();
        List<String> lines = FileReader.build().readLines(new File(hospitalAddressCsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_SEPARATOR, HOSPITAL_ADDRESS_FIELD_COUNT);
            if (parts.length == HOSPITAL_ADDRESS_FIELD_COUNT) {
                HospitalAddress hospital = ImmutableHospitalAddress.of(parts[HOSPITAL_ADDRESS_ID_COLUMN],
                        parts[HOSPITAL_ADDRESS_NAME_COLUMN],
                        parts[HOSPITAL_ADDRESS_ZIP_COLUMN],
                        parts[HOSPITAL_ADDRESS_CITY_COLUMN]);

                hospitalAddress.put(parts[SAMPLE_MAPPING_ID_COLUMN], hospital);
            } else {
                LOGGER.warn("Could not properly parse line in hospital adress csv: '{}'", line);
            }
        }
        return hospitalAddress;
    }
}
