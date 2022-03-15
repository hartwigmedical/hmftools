package com.hartwig.hmftools.patientreporter.panel;

import java.util.List;

import com.hartwig.hmftools.patientreporter.QsFormNumber;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailType;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum PanelFailReason {
    PANEL_FAILURE("panel_failure", false, QsFormNumber.FOR_102);


    @NotNull
    private final String identifier;
    private final boolean deepWGSDataAvailable;
    @NotNull
    private final QsFormNumber qsFormNumber;

    PanelFailReason(@NotNull final String identifier, final boolean deepWGSDataAvailable,
            @NotNull QsFormNumber qsFormNumber) {
        this.identifier = identifier;
        this.deepWGSDataAvailable = deepWGSDataAvailable;
        this.qsFormNumber = qsFormNumber;
    }

    @NotNull
    public String identifier() {
        return identifier;
    }

    public boolean isDeepWGSDataAvailable() {
        return deepWGSDataAvailable;
    }

    @NotNull
    public String qcFormNumber() {
        return qsFormNumber.display();
    }

    @Nullable
    public static PanelFailReason fromIdentifier(@Nullable String identifier) {
        if (identifier == null) {
            return null;
        }

        for (PanelFailReason reason : PanelFailReason.values()) {
            if (reason.identifier().equals(identifier)) {
                return reason;
            }
        }

        return null;
    }

    @NotNull
    public static List<String> validIdentifiers() {
        List<String> identifiers = Lists.newArrayList();
        for (PanelFailReason reason : PanelFailReason.values()) {
            identifiers.add(reason.identifier);
        }
        return identifiers;
    }
}
