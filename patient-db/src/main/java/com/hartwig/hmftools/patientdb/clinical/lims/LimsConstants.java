package com.hartwig.hmftools.patientdb.clinical.lims;

import java.time.format.DateTimeFormatter;

final class LimsConstants {

    static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    static final double DNA_MICRO_LITERS = 50D;

    private LimsConstants() {
    }
}
