package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CLINICALFINDINGS;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.ValidationFinding;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

class ValidationFindingDAO {

    @NotNull
    private final DSLContext context;

    ValidationFindingDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void clear() {
        context.truncate(CLINICALFINDINGS).execute();
    }

    void write(@NotNull List<ValidationFinding> findings) {
        context.batch(findings.stream()
                .map(finding -> context.insertInto(CLINICALFINDINGS,
                        CLINICALFINDINGS.LEVEL,
                        CLINICALFINDINGS.PATIENTID,
                        CLINICALFINDINGS.FORMSTATUS,
                        CLINICALFINDINGS.FORMLOCKED,
                        CLINICALFINDINGS.MESSAGE,
                        CLINICALFINDINGS.DETAILS)
                        .values(finding.level(),
                                finding.patientIdentifier(),
                                finding.formStatus().stateString(),
                                finding.formStatus().lockedString(),
                                finding.message(),
                                finding.details()))
                .collect(Collectors.toList())).execute();
    }
}
