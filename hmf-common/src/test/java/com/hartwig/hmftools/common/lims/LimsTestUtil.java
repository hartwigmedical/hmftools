package com.hartwig.hmftools.common.lims;

import java.time.LocalDate;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

final class LimsTestUtil {

    private LimsTestUtil() {
    }

    @NotNull
    static LocalDate toDate(@NotNull final String date) {
        return LocalDate.parse(date, LimsConstants.DATE_FORMATTER);
    }

    @NotNull
    static ImmutableLimsJsonSampleData.Builder createLimsSampleDataBuilder() {
        return ImmutableLimsJsonSampleData.builder()
                .sampleId(Strings.EMPTY)
                .tumorBarcodeId(Strings.EMPTY)
                .refBarcodeId(Strings.EMPTY)
                .arrivalDateString(Strings.EMPTY)
                .dnaConcentration(Strings.EMPTY)
                .primaryTumor(Strings.EMPTY)
                .labSopVersions(Strings.EMPTY)
                .labelSample(Strings.EMPTY)
                .projectName(Strings.EMPTY)
                .submission(Strings.EMPTY)
                .patientId(Strings.EMPTY)
                .requesterName(Strings.EMPTY)
                .requesterEmail(Strings.EMPTY)
                .shallowSeq(0);
    }
}
