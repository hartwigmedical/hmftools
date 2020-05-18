package com.hartwig.hmftools.common.lims.hospital;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HospitalModelFactory {

    private static final Logger LOGGER = LogManager.getLogger(HospitalModelFactory.class);

    private static final String HOSPITALS_ADDRESS_TSV = "hospital_address.tsv";
    private static final String HOSPITAL_CPCT_TSV = "hospital_cpct.tsv";
    private static final String HOSPITAL_DRUP_TSV = "hospital_drup.tsv";
    private static final String HOSPITAL_WIDE_TSV = "hospital_wide.tsv";
    private static final String SAMPLE_HOSPITAL_MAPPING_TSV = "sample_hospital_mapping.tsv";

    private static final int HOSPITAL_ADDRESS_ID_COLUMN = 0;
    private static final int HOSPITAL_ADDRESS_NAME_COLUMN = 1;
    private static final int HOSPITAL_ADDRESS_ZIP_COLUMN = 2;
    private static final int HOSPITAL_ADDRESS_CITY_COLUMN = 3;
    private static final int HOSPITAL_ADDRESS_FIELD_COUNT = 4;

    private static final int HOSPITAL_CONTACT_ID_COLUMN = 0;
    private static final int HOSPITAL_CONTACT_PI_COLUMN = 1;
    private static final int HOSPITAL_CONTACT_REQUESTER_NAME_COLUMN = 2;
    private static final int HOSPITAL_CONTACT_REQUESTER_EMAIL_COLUMN = 3;
    private static final int HOSPITAL_CONTACT_FIELD_COUNT_WIDE = 4;
    private static final int HOSPITAL_CONTACT_FIELD_COUNT_CPCT_DRUP = 2;

    private static final int SAMPLE_MAPPING_ID_COLUMN = 0;
    private static final int HOSPITAL_MAPPING_COLUMN = 1;
    private static final int FIELD_COUNT_SAMPLE_HOSPITAL_MAPPING = 2;

    private static final String FIELD_SEPARATOR = "\t";

    private HospitalModelFactory() {
    }

    @NotNull
    public static HospitalModel fromLimsDirectory(@NotNull String limsDirectory) throws IOException {
        String hospitalAddressTsv = limsDirectory + File.separator + HOSPITALS_ADDRESS_TSV;

        String hospitalContactCPCTTsv = limsDirectory + File.separator + HOSPITAL_CPCT_TSV;
        String hospitalContactDRUPTsv = limsDirectory + File.separator + HOSPITAL_DRUP_TSV;
        String hospitalContactWIDETsv = limsDirectory + File.separator + HOSPITAL_WIDE_TSV;

        String sampleHospitalMappingTsv = limsDirectory + File.separator + SAMPLE_HOSPITAL_MAPPING_TSV;

        Map<String, HospitalAddress> hospitalAddress = readFromHospitalAddress(hospitalAddressTsv);
        Map<String, HospitalPersons> hospitalContactCPCT =
                readFromHospitalContact(hospitalContactCPCTTsv, HOSPITAL_CONTACT_FIELD_COUNT_CPCT_DRUP);
        Map<String, HospitalPersons> hospitalContactDRUP =
                readFromHospitalContact(hospitalContactDRUPTsv, HOSPITAL_CONTACT_FIELD_COUNT_CPCT_DRUP);
        Map<String, HospitalPersons> hospitalContactWIDE =
                readFromHospitalContact(hospitalContactWIDETsv, HOSPITAL_CONTACT_FIELD_COUNT_WIDE);
        Map<String, String> sampleHospitalMapping = readFromSampleToHospitalMapping(sampleHospitalMappingTsv);

        HospitalModel hospitalModel = ImmutableHospitalModel.builder()
                .hospitalAddress(hospitalAddress)
                .sampleToHospitalMapping(sampleHospitalMapping)
                .hospitalContactCPCT(hospitalContactCPCT)
                .hospitalContactDRUP(hospitalContactDRUP)
                .hospitalContactWIDE(hospitalContactWIDE)
                .build();

        validateModelIntegrity(hospitalModel);

        return hospitalModel;
    }

    @VisibleForTesting
    static boolean validateModelIntegrity(@NotNull HospitalModel hospitalModel) {
        Set<String> hospitalIdsInAddressList = hospitalModel.hospitalAddress().keySet();
        Set<String> keyCPCT = hospitalModel.hospitalContactCPCT().keySet();
        Set<String> keyDRUP = hospitalModel.hospitalContactDRUP().keySet();
        Set<String> keyWIDE = hospitalModel.hospitalContactWIDE().keySet();
        Set<String> keySampleMapping = Sets.newHashSet(hospitalModel.sampleToHospitalMapping().values());

        boolean allCorrect = true;
        for (String CPCT : keyCPCT) {
            if (!hospitalIdsInAddressList.contains(CPCT)) {
                allCorrect = false;
                LOGGER.warn("CPCT hospital ID is not present in hospital address list: '{}'", CPCT);
            }
        }

        for (String DRUP : keyDRUP) {
            if (!hospitalIdsInAddressList.contains(DRUP)) {
                allCorrect = false;
                LOGGER.warn("DRUP hospital ID is not present in hospital address list: '{}'", DRUP);
            }
        }

        for (String WIDE : keyWIDE) {
            if (!hospitalIdsInAddressList.contains(WIDE)) {
                allCorrect = false;
                LOGGER.warn("WIDE hospital ID is not present in hospital address list: '{}'", WIDE);
            }
        }

        for (String sampleMapping : keySampleMapping) {
            if (!hospitalIdsInAddressList.contains(sampleMapping)) {
                allCorrect = false;
                LOGGER.warn("Sample mapping hospital ID is not present in hospital address list: '{}'", sampleMapping);
            }
        }

        return allCorrect;
    }

    @NotNull
    @VisibleForTesting
    static Map<String, HospitalAddress> readFromHospitalAddress(@NotNull String hospitalAddressTsv) throws IOException {
        Map<String, HospitalAddress> hospitalAddressMap = Maps.newHashMap();
        List<String> lines = FileReader.build().readLines(new File(hospitalAddressTsv).toPath());
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
    static Map<String, HospitalPersons> readFromHospitalContact(@NotNull String hospitalContactTsv, int expectedFieldCount)
            throws IOException {
        Map<String, HospitalPersons> hospitalContactMap = Maps.newHashMap();
        List<String> lines = FileReader.build().readLines(new File(hospitalContactTsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_SEPARATOR);
            if (parts.length == expectedFieldCount) {
                String requesterName = parts.length > 2 ? parts[HOSPITAL_CONTACT_REQUESTER_NAME_COLUMN] : null;
                String requesterEmail = parts.length > 2 ? parts[HOSPITAL_CONTACT_REQUESTER_EMAIL_COLUMN] : null;

                HospitalPersons hospitalPersons = ImmutableHospitalPersons.builder()
                        .hospitalPI(parts[HOSPITAL_CONTACT_PI_COLUMN])
                        .requesterName(requesterName)
                        .requesterEmail(requesterEmail)
                        .build();

                hospitalContactMap.put(parts[HOSPITAL_CONTACT_ID_COLUMN], hospitalPersons);
            } else {
                LOGGER.warn("Could not properly parse line in hospital contact tsv: '{}'", line);
            }
        }
        return hospitalContactMap;
    }

    @NotNull
    @VisibleForTesting
    static Map<String, String> readFromSampleToHospitalMapping(@NotNull String sampleHospitalMappingTsv) throws IOException {
        Map<String, String> hospitalPerSampleMap = Maps.newHashMap();

        List<String> lines = FileReader.build().readLines(new File(sampleHospitalMappingTsv).toPath());
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
