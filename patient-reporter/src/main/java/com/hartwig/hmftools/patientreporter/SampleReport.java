package com.hartwig.hmftools.patientreporter;

import java.time.LocalDate;

import com.hartwig.hmftools.common.clinical.PatientTumorLocation;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.lims.hospital.HospitalContactData;
import com.hartwig.hmftools.patientreporter.cfreport.data.DataUtil;

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
    public abstract PatientTumorLocation patientTumorLocation();

    @NotNull
    public abstract LimsGermlineReportingLevel germlineReportingLevel();

    public abstract boolean reportViralInsertions();

    @Nullable
    public abstract LocalDate refArrivalDate();

    @Nullable
    public abstract LocalDate tumorArrivalDate();

    @NotNull
    public abstract String shallowSeqPurityString();

    @NotNull
    public abstract String labProcedures();

    @NotNull
    public abstract String cohort();

    @NotNull
    public abstract String projectName();

    @NotNull
    public abstract String submissionId();

    @NotNull
    public abstract HospitalContactData hospitalContactData();

    @NotNull
    public abstract String hospitalPatientId();

    @Nullable
    public abstract String hospitalPathologySampleId();

    @Nullable
    @Value.Derived
    public String refSampleBarcode() {
        return sampleMetadata().refSampleBarcode();
    }

    @NotNull
    @Value.Derived
    public String tumorSampleId() {
        return sampleMetadata().tumorSampleId();
    }

    @NotNull
    @Value.Derived
    public String tumorSampleBarcode() {
        return sampleMetadata().tumorSampleBarcode();
    }

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

    @NotNull
    @Value.Derived
    public String primaryTumorLocationString() {
        PatientTumorLocation entry = patientTumorLocation();
        if (entry != null) {
            if (!entry.primaryTumorSubLocation().isEmpty()) {
                return entry.primaryTumorLocation() + " (" + entry.primaryTumorSubLocation() + ")";
            } else {
                return entry.primaryTumorLocation();
            }
        } else {
            return Strings.EMPTY;
        }
    }

    @NotNull
    @Value.Derived
    public String primaryTumorTypeString() {
        PatientTumorLocation entry = patientTumorLocation();
        if (entry != null) {
            if (!entry.primaryTumorSubType().isEmpty()) {
                return entry.primaryTumorSubType();
            } else {
                return entry.primaryTumorType();
            }
        } else {
            return Strings.EMPTY;
        }
    }

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
