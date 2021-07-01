package com.hartwig.hmftools.common.lims.hospital;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HospitalModelFactory {

    private static final Logger LOGGER = LogManager.getLogger(HospitalModelFactory.class);

    private static final String HOSPITALS_ADDRESS_TSV = "hospital_address.tsv";
    private static final String HOSPITAL_CPCT_TSV = "hospital_cpct.tsv";
    private static final String HOSPITAL_DRUP_TSV = "hospital_drup.tsv";
    private static final String HOSPITAL_WIDE_TSV = "hospital_wide.tsv";
    private static final String HOSPITAL_COREDB_TSV = "hospital_coredb.tsv";
    private static final String HOSPITAL_ACTIN_TSV = "hospital_actin.tsv";
    private static final String SAMPLE_HOSPITAL_MAPPING_TSV = "sample_hospital_mapping.tsv";

    private static final int HOSPITAL_ADDRESS_ID_COLUMN = 0;
    private static final int HOSPITAL_ADDRESS_NAME_COLUMN = 1;
    private static final int HOSPITAL_ADDRESS_ZIP_COLUMN = 2;
    private static final int HOSPITAL_ADDRESS_CITY_COLUMN = 3;
    private static final int HOSPITAL_ADDRESS_FIELD_COUNT = 4;

    private static final int HOSPITAL_PERSONS_ID_COLUMN = 0;
    private static final int HOSPITAL_PERSONS_PI_COLUMN = 1;
    private static final int HOSPITAL_PERSONS_REQUESTER_NAME_COLUMN = 2;
    private static final int HOSPITAL_PERSONS_REQUESTER_EMAIL_COLUMN = 3;
    private static final int HOSPITAL_PERSONS_FIELD_COUNT_WIDE = 4;
    private static final int HOSPITAL_PERSONS_FIELD_COUNT_COREDB = 4;
    private static final int HOSPITAL_PERSONS_FIELD_COUNT_CPCT_DRUP_ACTIN = 2;

    private static final int SAMPLE_MAPPING_ID_COLUMN = 0;
    private static final int HOSPITAL_MAPPING_COLUMN = 1;
    private static final int FIELD_COUNT_SAMPLE_HOSPITAL_MAPPING = 2;

    private static final String FIELD_SEPARATOR = "\t";

    private HospitalModelFactory() {
    }

    @NotNull
    public static HospitalModel fromLimsDirectory(@NotNull String limsDirectory) throws IOException {
        String hospitalAddressTsv = limsDirectory + File.separator + HOSPITALS_ADDRESS_TSV;

        String hospitalPersonsCPCTTsv = limsDirectory + File.separator + HOSPITAL_CPCT_TSV;
        String hospitalPersonsDRUPTsv = limsDirectory + File.separator + HOSPITAL_DRUP_TSV;
        String hospitalPersonsWIDETsv = limsDirectory + File.separator + HOSPITAL_WIDE_TSV;
        String hospitalPersonsCOREDBTsv = limsDirectory + File.separator + HOSPITAL_COREDB_TSV;

        String sampleHospitalMappingTsv = limsDirectory + File.separator + SAMPLE_HOSPITAL_MAPPING_TSV;

        Map<String, HospitalAddress> hospitalAddressMap = readFromHospitalAddressTsv(hospitalAddressTsv);
        Map<String, HospitalPersons> hospitalPersonsCPCT =
                readFromHospitalPersonsTsv(hospitalPersonsCPCTTsv, HOSPITAL_PERSONS_FIELD_COUNT_CPCT_DRUP_ACTIN, "CPCT");
        Map<String, HospitalPersons> hospitalPersonsDRUP =
                readFromHospitalPersonsTsv(hospitalPersonsDRUPTsv, HOSPITAL_PERSONS_FIELD_COUNT_CPCT_DRUP_ACTIN, "DRUP");
        Map<String, HospitalPersons> hospitalPersonsWIDE =
                readFromHospitalPersonsTsv(hospitalPersonsWIDETsv, HOSPITAL_PERSONS_FIELD_COUNT_WIDE, "WIDE");
        Map<String, HospitalPersons> hospitalPersonsCOREDB =
                readFromHospitalPersonsTsv(hospitalPersonsCOREDBTsv, HOSPITAL_PERSONS_FIELD_COUNT_COREDB, "COREDB");
        Map<String, HospitalPersons> hospitalPersonsACTIN =
                readFromHospitalPersonsTsv(hospitalPersonsCOREDBTsv, HOSPITAL_PERSONS_FIELD_COUNT_CPCT_DRUP_ACTIN, "ACTIN");
        Map<String, String> sampleHospitalMapping = readFromSampleToHospitalMappingTsv(sampleHospitalMappingTsv);

        HospitalModel hospitalModel = ImmutableHospitalModel.builder()
                .hospitalAddressMap(hospitalAddressMap)
                .hospitalPersonsCPCT(hospitalPersonsCPCT)
                .hospitalPersonsDRUP(hospitalPersonsDRUP)
                .hospitalPersonsWIDE(hospitalPersonsWIDE)
                .hospitalPersonsCOREDB(hospitalPersonsCOREDB)
                .hospitalPersonsACTIN(hospitalPersonsACTIN)
                .sampleToHospitalMapping(sampleHospitalMapping)
                .build();

        HospitalChecker.validateModelIntegrity(hospitalModel);

        return hospitalModel;
    }

    @NotNull
    @VisibleForTesting
    static Map<String, HospitalAddress> readFromHospitalAddressTsv(@NotNull String hospitalAddressTsv) throws IOException {
        Map<String, HospitalAddress> hospitalAddressMap = Maps.newHashMap();
        List<String> lines = Files.readAllLines(new File(hospitalAddressTsv).toPath());

        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_SEPARATOR);
            if (parts.length == HOSPITAL_ADDRESS_FIELD_COUNT) {
                HospitalAddress hospitalAddress = ImmutableHospitalAddress.builder()
                        .hospitalName(parts[HOSPITAL_ADDRESS_NAME_COLUMN])
                        .hospitalZip(parts[HOSPITAL_ADDRESS_ZIP_COLUMN])
                        .hospitalCity(parts[HOSPITAL_ADDRESS_CITY_COLUMN])
                        .build();

                hospitalAddressMap.put(parts[HOSPITAL_ADDRESS_ID_COLUMN], hospitalAddress);
            } else {
                LOGGER.warn("Could not properly parse line in hospital address tsv: '{}'", line);
            }
        }
        return hospitalAddressMap;
    }

    @NotNull
    @VisibleForTesting
    static Map<String, HospitalPersons> readFromHospitalPersonsTsv(@NotNull String hospitalPersonsTsv, int expectedFieldCount,
            @NotNull String cohort) throws IOException {
        Map<String, HospitalPersons> hospitalPersonsMap = Maps.newHashMap();
        List<String> lines = Files.readAllLines(new File(hospitalPersonsTsv).toPath());

        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_SEPARATOR);
            if (parts.length == expectedFieldCount) {
                String requesterName = parts.length > 2 ? parts[HOSPITAL_PERSONS_REQUESTER_NAME_COLUMN] : null;
                String requesterEmail = parts.length > 2 ? parts[HOSPITAL_PERSONS_REQUESTER_EMAIL_COLUMN] : null;

                HospitalPersons hospitalPersons = ImmutableHospitalPersons.builder()
                        .hospitalPI(cohort.equals("COREDB") ? requesterName : parts[HOSPITAL_PERSONS_PI_COLUMN])
                        .requesterName(requesterName)
                        .requesterEmail(requesterEmail)
                        .build();

                hospitalPersonsMap.put(parts[HOSPITAL_PERSONS_ID_COLUMN], hospitalPersons);
            } else {
                LOGGER.warn("Could not properly parse line in hospital persons tsv: '{}'", line);
            }
        }
        return hospitalPersonsMap;
    }

    @NotNull
    @VisibleForTesting
    static Map<String, String> readFromSampleToHospitalMappingTsv(@NotNull String sampleHospitalMappingTsv) throws IOException {
        Map<String, String> hospitalPerSampleMap = Maps.newHashMap();
        List<String> lines = Files.readAllLines(new File(sampleHospitalMappingTsv).toPath());

        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_SEPARATOR);
            if (parts.length == FIELD_COUNT_SAMPLE_HOSPITAL_MAPPING) {
                hospitalPerSampleMap.put(parts[SAMPLE_MAPPING_ID_COLUMN], parts[HOSPITAL_MAPPING_COLUMN]);
            } else {
                LOGGER.warn("Could not properly parse line in sample hospital mapping tsv: '{}'", line);
            }
        }
        return hospitalPerSampleMap;
    }
}
