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

    private static final String HOSPITALS_CSV = "hospitals.csv";
    private static final String SAMPLE_HOSPITAL_MAPPING_CSV = "sample_hospital_mapping.csv";
    private static final String HOSPITALS_CORE_CSV = "hospitals_core.csv";

    private static final int HOSPITAL_ID_COLUMN = 0;
    private static final int INTERNAL_HOSPITAL_NAME_COLUMN = 1;
    private static final int CPCT_PI_COLUMN = 3;
    private static final int CPCT_RECIPIENTS_COLUMN = 4;
    private static final int DRUP_PI_COLUMN = 6;
    private static final int DRUP_RECIPIENTS_COLUMN = 7;
    private static final int WIDE_PI_COLUMN = 9;
    private static final int WIDE_RECIPIENTS_COLUMN = 10;
    private static final int EXTERNAL_HOSPITAL_NAME_COLUMN = 13;
    private static final int ADDRESS_ZIP_COLUMN = 14;
    private static final int ADDRESS_CITY_COLUMN = 15;
    private static final int HOSPITAL_FIELD_COUNT = 16;

    private static final int SAMPLE_MAPPING_ID_COLUMN = 0;
    private static final int HOSPITAL_MAPPING_COLUMN = 1;
    private static final int FIELD_COUNT_SAMPLE_HOSPITAL_MAPPING = 2;

    private static final int HOSPITAL_CORE_ID_COLUMN = 0;
    private static final int INTERNAL_CORE_HOSPITAL_NAME_COLUMN = 1;
    private static final int EXTERNAL_CORE_HOSPITAL_NAME_COLUMN = 2;
    private static final int ADDRESS_CORE_ZIP_COLUMN = 3;
    private static final int ADDRESS_CORE_CITY_COLUMN = 4;
    private static final int HOSPITAL_CORE_FIELD_COUNT = 5;

    private static final String FIELD_SEPARATOR = ",";

    private HospitalModelFactory() {
    }

    @NotNull
    public static HospitalModel fromHospitalDirectory(@NotNull String hospitalDirectory) throws IOException {
        String hospitalCsvPath = hospitalDirectory + File.separator + HOSPITALS_CSV;
        String sampleHospitalMappingCsv = hospitalDirectory + File.separator + SAMPLE_HOSPITAL_MAPPING_CSV;
        String hospitalCoreCsv = hospitalDirectory + File.separator + HOSPITALS_CORE_CSV;

        Map<String, HospitalData> hospitalPerId = readFromHospitalCsv(hospitalCsvPath);
        Map<String, HospitalSampleMapping> sampleHospitalMapping = readFromSampleHospitalMapping(sampleHospitalMappingCsv);
        Map<String, HospitalCore> hospitalCoreMap = readFromCoreHospitalMapping(hospitalCoreCsv);

        return ImmutableHospitalModel.of(hospitalPerId, sampleHospitalMapping, hospitalCoreMap, Maps.newHashMap(), Maps.newHashMap());
    }

    @NotNull
    public static HospitalModel empty() {
        return ImmutableHospitalModel.of(Maps.newHashMap(), Maps.newHashMap(), Maps.newHashMap(), Maps.newHashMap(), Maps.newHashMap());
    }

    @NotNull
    private static Map<String, HospitalCore> readFromCoreHospitalMapping(@NotNull String hospitalCsv) throws IOException {
        Map<String, HospitalCore> hospitalCore = Maps.newHashMap();

        List<String> lines = FileReader.build().readLines(new File(hospitalCsv).toPath());
        for (String line : lines) {
            String[] parts = line.split(FIELD_SEPARATOR, HOSPITAL_CORE_FIELD_COUNT);
            if (parts.length == HOSPITAL_CORE_FIELD_COUNT) {
                HospitalCore hospital = ImmutableHospitalCore.of(parts[INTERNAL_CORE_HOSPITAL_NAME_COLUMN],
                        parts[EXTERNAL_CORE_HOSPITAL_NAME_COLUMN],
                        parts[ADDRESS_CORE_ZIP_COLUMN],
                        parts[ADDRESS_CORE_CITY_COLUMN]);

                hospitalCore.put(parts[HOSPITAL_CORE_ID_COLUMN], hospital);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in hospital csv: '{}'", line);
            }
        }

        return hospitalCore;
    }

    @NotNull
    private static Map<String, HospitalData> readFromHospitalCsv(@NotNull String hospitalCsv) throws IOException {
        Map<String, HospitalData> hospitalPerId = Maps.newHashMap();

        List<String> lines = FileReader.build().readLines(new File(hospitalCsv).toPath());
        for (String line : lines) {
            String[] parts = line.split(FIELD_SEPARATOR, HOSPITAL_FIELD_COUNT);
            if (parts.length == HOSPITAL_FIELD_COUNT) {
                HospitalData hospital = ImmutableHospitalData.of(parts[INTERNAL_HOSPITAL_NAME_COLUMN],
                        parts[CPCT_RECIPIENTS_COLUMN],
                        parts[DRUP_RECIPIENTS_COLUMN],
                        parts[WIDE_RECIPIENTS_COLUMN],
                        parts[EXTERNAL_HOSPITAL_NAME_COLUMN],
                        parts[ADDRESS_ZIP_COLUMN],
                        parts[ADDRESS_CITY_COLUMN],
                        parts[CPCT_PI_COLUMN],
                        parts[DRUP_PI_COLUMN],
                        parts[WIDE_PI_COLUMN]);

                hospitalPerId.put(parts[HOSPITAL_ID_COLUMN], hospital);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in hospital csv: '{}'", line);
            }
        }

        return hospitalPerId;
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
}
