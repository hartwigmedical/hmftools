package com.hartwig.hmftools.patientreporter;

import java.util.List;
import java.util.Optional;

import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlCData;
import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlProperty;
import com.hartwig.hmftools.common.peach.PeachGenotype;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public interface PatientReport {

    @NotNull
    SampleReport sampleReport();

    @JacksonXmlProperty(isAttribute = true, localName = "user")
    @NotNull
    default String user() {
        String systemUser = System.getProperty("user.name");
        String userName = Strings.EMPTY;
        String trainedEmployee = " (trained IT employee)";
        String combinedUserName = Strings.EMPTY;
        if (systemUser.equals("lieke") || systemUser.equals("liekeschoenmaker") || systemUser.equals("lschoenmaker")) {
            userName = "Lieke Schoenmaker";
            combinedUserName = userName + trainedEmployee;
        } else if (systemUser.equals("korneel") || systemUser.equals("korneelduyvesteyn") || systemUser.equals("kduyvesteyn")) {
            userName = "Korneel Duyvesteyn";
            combinedUserName = userName + trainedEmployee;
        } else if (systemUser.equals("sandra") || systemUser.equals("sandravandenbroek") || systemUser.equals("sandravdbroek")
                || systemUser.equals("s_vandenbroek") || systemUser.equals("svandenbroek")) {
            userName = "Sandra van den Broek";
            combinedUserName = userName + trainedEmployee;
        } else if (systemUser.equals("daphne") || systemUser.equals("d_vanbeek") || systemUser.equals("daphnevanbeek")
                || systemUser.equals("dvanbeek")) {
            userName = "Daphne van Beek";
            combinedUserName = userName + trainedEmployee;
        } else if (systemUser.equals("root")) {
            combinedUserName = "automatically";
        } else {
            userName = systemUser;
            combinedUserName = userName + trainedEmployee;
        }

        if (combinedUserName.endsWith(trainedEmployee)) {
            combinedUserName = "by " + combinedUserName;
        }

        return combinedUserName + " and checked by a trained Clinical Molecular Biologist in Pathology (KMBP)";

    }

    @JacksonXmlProperty(isAttribute = true, localName = "qsFormNumber")
    @NotNull
    String qsFormNumber();

    @JacksonXmlProperty(isAttribute = true, localName = "comments")
    @NotNull
    Optional<String> comments();

    @JacksonXmlProperty(isAttribute = true, localName = "isCorrectedReport")
    boolean isCorrectedReport();

    @JacksonXmlProperty(isAttribute = true, localName = "isCorrectedReportExtern")
    boolean isCorrectedReportExtern();

    @JacksonXmlProperty(isAttribute = true, localName = "signaturePath")
    @NotNull
    String signaturePath();

    @JacksonXmlProperty(isAttribute = true, localName = "logoRVAPath")
    @NotNull
    String logoRVAPath();

    @JacksonXmlProperty(isAttribute = true, localName = "logoCompanyPath")
    @NotNull
    String logoCompanyPath();

    @JacksonXmlProperty(isAttribute = true, localName = "udiDi")
    @NotNull
    String udiDi();

    @JacksonXmlProperty(localName = "peachGenotypes")
    @JacksonXmlCData
    @NotNull
    List<PeachGenotype> peachGenotypes();

    @JacksonXmlProperty(isAttribute = true, localName = "reportDate")
    @NotNull
    String reportDate();

    @JacksonXmlProperty(isAttribute = true, localName = "isWGSreport")
    boolean isWGSreport();
}
