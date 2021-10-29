package com.hartwig.hmftools.common.lims.hospital;

import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HospitalChecker {

    private static final Logger LOGGER = LogManager.getLogger(HospitalChecker.class);

    private HospitalChecker() {
    }

    @VisibleForTesting
    static boolean validateModelIntegrity(@NotNull HospitalModel hospitalModel) {
        Set<String> hospitalIdsInAddressList = hospitalModel.hospitalAddressMap().keySet();
        Set<String> keyCPCT = hospitalModel.hospitalPersonsCPCT().keySet();
        Set<String> keyDRUP = hospitalModel.hospitalPersonsDRUP().keySet();
        Set<String> keyWIDE = hospitalModel.hospitalPersonsWIDE().keySet();
        Set<String> keyCOREDB = hospitalModel.hospitalPersonsCOREDB().keySet();
        Set<String> keyACTIN = hospitalModel.hospitalPersonsACTIN().keySet();
        Set<String> keyGLOW = hospitalModel.hospitalPersonsGLOW().keySet();
        Set<String> keyOPTIC = hospitalModel.hospitalPersonsOPTIC().keySet();
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

        for (String COREDB : keyCOREDB) {
            if (!hospitalIdsInAddressList.contains(COREDB)) {
                allCorrect = false;
                LOGGER.warn("COREDB hospital ID is not present in hospital address list: '{}'", COREDB);
            }
        }

        for (String ACTIN : keyACTIN) {
            if (!hospitalIdsInAddressList.contains(ACTIN)) {
                allCorrect = false;
                LOGGER.warn("ACTIN hospital ID is not present in hospital address list: '{}'", ACTIN);
            }
        }

        for (String GLOW : keyGLOW) {
            if (!hospitalIdsInAddressList.contains(GLOW)) {
                allCorrect = false;
                LOGGER.warn("ACTIN hospital ID is not present in hospital address list: '{}'", GLOW);
            }
        }

        for (String OPTIC : keyOPTIC) {
            if (!hospitalIdsInAddressList.contains(OPTIC)) {
                allCorrect = false;
                LOGGER.warn("ACTIN hospital ID is not present in hospital address list: '{}'", OPTIC);
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
}
