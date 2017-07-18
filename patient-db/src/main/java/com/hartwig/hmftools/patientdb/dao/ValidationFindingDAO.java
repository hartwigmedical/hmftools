package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CLINICALLOGS;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

public class ValidationFindingDAO {
    @NotNull
    private final DSLContext context;

    ValidationFindingDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void clear() {
        context.truncate(CLINICALLOGS).execute();
    }

    void write(@NotNull final List<ValidationFinding> findings) {
        context.batch(findings.stream()
                .map(finding -> context.insertInto(CLINICALLOGS, CLINICALLOGS.PATIENTID, CLINICALLOGS.ECRFITEM, CLINICALLOGS.FORMSTATUS,
                        CLINICALLOGS.FORMLOCKED, CLINICALLOGS.MESSAGE)
                        .values(finding.patientId(), finding.ecrfItem(), finding.formStatus(), finding.formLocked(), finding.message()))
                .collect(Collectors.toList())).execute();
    }
}
