package com.hartwig.hmftools.patientreporter;

import java.time.LocalDate;

import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlProperty;
import com.fasterxml.jackson.dataformat.xml.annotation.JacksonXmlRootElement;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.lims.hospital.HospitalContactData;
import com.hartwig.hmftools.common.utils.DataUtil;

import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class SampleReport {

    @NotNull
    public abstract SampleMetadata sampleMetadata();

    @Nullable
    public abstract PatientPrimaryTumor patientPrimaryTumor();

    @JacksonXmlProperty(isAttribute = true, localName = "biopsyLocation")
    @Nullable
    public abstract String biopsyLocation();

    @JacksonXmlProperty(isAttribute = true, localName = "germlineReportingLevel")
    @NotNull
    public abstract LimsGermlineReportingLevel germlineReportingLevel();

    @JacksonXmlProperty(isAttribute = true, localName = "reportViralPresence")
    public abstract boolean reportViralPresence();

    @JacksonXmlProperty(isAttribute = true, localName = "reportPharmogenetics")
    public abstract boolean reportPharmogenetics();

    @JacksonXmlProperty(isAttribute = true, localName = "refArrivalDate")
    @Nullable
    public abstract LocalDate refArrivalDate();

    @JacksonXmlProperty(isAttribute = true, localName = "tumorArrivalDate")
    @Nullable
    public abstract LocalDate tumorArrivalDate();

    @JacksonXmlProperty(isAttribute = true, localName = "shallowSeqPurityString")
    @NotNull
    public abstract String shallowSeqPurityString();

    @JacksonXmlProperty(isAttribute = true, localName = "labProcedures")
    @NotNull
    public abstract String labProcedures();

    @NotNull
    public abstract LimsCohortConfig cohort();

    @JacksonXmlProperty(isAttribute = true, localName = "projectName")
    @NotNull
    public abstract String projectName();

    @JacksonXmlProperty(isAttribute = true, localName = "submissionId")
    @NotNull
    public abstract String submissionId();

    @NotNull
    public abstract HospitalContactData hospitalContactData();

    @JacksonXmlProperty(isAttribute = true, localName = "hospitalPatientId")
    @NotNull
    public abstract String hospitalPatientId();

    @JacksonXmlProperty(isAttribute = true, localName = "hospitalPathologySampleId")
    @Nullable
    public abstract String hospitalPathologySampleId();

    @JacksonXmlProperty(isAttribute = true, localName = "refSampleBarcode")
    @Nullable
    @Value.Derived
    public String refSampleBarcode() {
        return sampleMetadata().refSampleBarcode();
    }

    @JacksonXmlProperty(isAttribute = true, localName = "tumorSampleId")
    @NotNull
    @Value.Derived
    public String tumorSampleId() {
        return sampleMetadata().tumorSampleId();
    }

    @JacksonXmlProperty(isAttribute = true, localName = "tumorSampleBarcode")
    @NotNull
    @Value.Derived
    public String tumorSampleBarcode() {
        return sampleMetadata().tumorSampleBarcode();
    }

    @JacksonXmlProperty(isAttribute = true, localName = "sampleNameForReport")
    @NotNull
    @Value.Derived
    public String sampleNameForReport() {
        return sampleMetadata().sampleNameForReport();
    }

    @JacksonXmlProperty(isAttribute = true, localName = "earliestArrivalDate")
    @Nullable
    @Value.Derived
    public String earliestArrivalDate() {
        LocalDate refDate = refArrivalDate();
        LocalDate sampleDate = tumorArrivalDate();

        if (sampleDate == null) {
            return null;
        } else if (refDate == null || sampleDate.isBefore(refDate)) {
            return DataUtil.formatDate(sampleDate);
        } else {
            return DataUtil.formatDate(refDate);
        }
    }

    @JacksonXmlProperty(isAttribute = true, localName = "biopsyLocationString")
    @NotNull
    @Value.Derived
    public String biopsyLocationString() {
        String biopsyLocation = biopsyLocation();
        if (biopsyLocation != null) {
            return biopsyLocation;
        } else {
            return Strings.EMPTY;
        }
    }

    @JacksonXmlProperty(isAttribute = true, localName = "primaryTumorLocationString")
    @NotNull
    @Value.Derived
    public String primaryTumorLocationString() {
        PatientPrimaryTumor entry = patientPrimaryTumor();
        if (entry != null) {
            if (!entry.subLocation().isEmpty()) {
                return entry.location() + " (" + entry.subLocation() + ")";
            } else {
                return entry.location();
            }
        } else {
            return Strings.EMPTY;
        }
    }

    @JacksonXmlProperty(isAttribute = true, localName = "primaryTumorTypeString")
    @NotNull
    @Value.Derived
    public String primaryTumorTypeString() {
        PatientPrimaryTumor entry = patientPrimaryTumor();
        if (entry != null) {
            if (!entry.subType().isEmpty()) {
                return entry.subType();
            } else {
                return entry.type();
            }
        } else {
            return Strings.EMPTY;
        }
    }

    @JacksonXmlProperty(isAttribute = true, localName = "addressee")
    @NotNull
    @Value.Derived
    public String addressee() {
        if (!hospitalContactData().hospitalPI().equals(Lims.NOT_AVAILABLE_STRING)) {
            return hospitalContactData().hospitalPI() + ", " + hospitalContactData().hospitalName() + ", "
                    + hospitalContactData().hospitalAddress();
        } else {
            return hospitalContactData().hospitalName() + ", " + hospitalContactData().hospitalAddress();
        }
    }
}
