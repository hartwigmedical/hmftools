package com.hartwig.hmftools.patientreporter;

import java.util.Optional;

import org.jetbrains.annotations.NotNull;

public interface PatientReport {

    @NotNull
    SampleReport sampleReport();

    @NotNull
    default String user() {
        String systemUser = System.getProperty("user.name");
        if (systemUser.equals("lieke") || systemUser.equals("liekeschoenmaker")) {
            return "Lieke Schoenmaker";
        } else if (systemUser.equals("korneel") || systemUser.equals("korneelduyvesteyn")) {
            return "Korneel Duyvesteyn";
        } else if (systemUser.equals("sandra") || systemUser.equals("sandravandenbroek")) {
            return "Sandra van den Broek";
        } else if (systemUser.equals("root")) {
            return "trained IT employee and trained clinical molecular biologist (KMBP)";
        } else {
            return systemUser;
        }
    }

    @NotNull
    String qsFormNumber();

    @NotNull
    Optional<String> comments();

    boolean isCorrectedReport();

    @NotNull
    String signaturePath();

    @NotNull
    String logoRVAPath();

    @NotNull
    String logoCompanyPath();
}
