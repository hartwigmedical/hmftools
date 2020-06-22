package com.hartwig.hmftools.patientreporter.qcfail;

import java.util.List;

import com.hartwig.hmftools.patientreporter.QsFormNumber;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum QCFailReason {
    TECHNICAL_FAILURE("technical_failure", QCFailType.TECHNICAL_FAILURE, false, QsFormNumber.FOR_102),
    SUFFICIENT_TCP_QC_FAILURE("sufficient_tcp_qc_failure", QCFailType.LOW_QUALITY_BIOPSY, true, QsFormNumber.FOR_083),
    INSUFFICIENT_TCP_SHALLOW_WGS("insufficient_tcp_shallow_wgs", QCFailType.LOW_QUALITY_BIOPSY, false, QsFormNumber.FOR_100),
    INSUFFICIENT_TCP_DEEP_WGS("insufficient_tcp_deep_wgs", QCFailType.LOW_QUALITY_BIOPSY, true, QsFormNumber.FOR_100),
    INSUFFICIENT_DNA("insufficient_dna", QCFailType.LOW_QUALITY_BIOPSY, false, QsFormNumber.FOR_082),
    UNDEFINED(Strings.EMPTY, QCFailType.UNDEFINED, false, QsFormNumber.FOR_UNDEFINED);

    @NotNull
    private final String identifier;
    @NotNull
    private final QCFailType type;
    private final boolean deepWGSDataAvailable;
    @NotNull
    private final QsFormNumber qsFormNumber;

    QCFailReason(@NotNull final String identifier, @NotNull final QCFailType type, final boolean deepWGSDataAvailable,
            @NotNull QsFormNumber qsFormNumber) {
        this.identifier = identifier;
        this.type = type;
        this.deepWGSDataAvailable = deepWGSDataAvailable;
        this.qsFormNumber = qsFormNumber;
    }

    @NotNull
    public String identifier() {
        return identifier;
    }

    @NotNull
    public QCFailType type() {
        return type;
    }

    public boolean isDeepWGSDataAvailable() {
        return deepWGSDataAvailable;
    }

    @NotNull
    public String qcFormNumber() {
        return qsFormNumber.display();
    }

    @NotNull
    public static QCFailReason fromIdentifier(@Nullable String identifier) {
        if (identifier == null) {
            return UNDEFINED;
        }

        for (QCFailReason reason : QCFailReason.values()) {
            if (reason.identifier().equals(identifier)) {
                return reason;
            }
        }

        return UNDEFINED;
    }

    @NotNull
    public static List<String> validIdentifiers() {
        List<String> identifiers = Lists.newArrayList();
        for (QCFailReason reason : QCFailReason.values()) {
            if (reason != QCFailReason.UNDEFINED) {
                identifiers.add(reason.identifier);
            }
        }
        return identifiers;
    }
}
