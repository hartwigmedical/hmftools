package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.HLATYPE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.HLATYPEDETAILS;

import java.time.LocalDateTime;
import java.util.List;

import com.hartwig.hmftools.common.hla.HlaType;
import com.hartwig.hmftools.common.hla.HlaTypeDetails;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep15;

class HlaTypeDAO {

    @NotNull
    private final DSLContext context;

    HlaTypeDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeType(@NotNull final String sample, @NotNull HlaType type) {
        context.delete(HLATYPE).where(HLATYPE.SAMPLEID.eq(sample)).execute();

        LocalDateTime timestamp = LocalDateTime.now();
        context.insertInto(HLATYPE,
                HLATYPE.MODIFIED,
                HLATYPE.SAMPLEID,
                HLATYPE.STATUS,
                HLATYPE.TYPEA1,
                HLATYPE.TYPEA2,
                HLATYPE.TYPEB1,
                HLATYPE.TYPEB2,
                HLATYPE.TYPEC1,
                HLATYPE.TYPEC2,
                HLATYPE.SOMATICVARIANTS)
                .values(timestamp,
                        sample,
                        type.status(),
                        type.typeA1(),
                        type.typeA2(),
                        type.typeB1(),
                        type.typeB2(),
                        type.typeC1(),
                        type.typeC2(),
                        type.somaticVariants())
                .execute();
    }

    void writeTypeDetails(@NotNull final String sample, @NotNull List<HlaTypeDetails> typeDetails) {
        context.delete(HLATYPEDETAILS).where(HLATYPEDETAILS.SAMPLEID.eq(sample)).execute();
        LocalDateTime timestamp = LocalDateTime.now();

        InsertValuesStep15 inserter = context.insertInto(HLATYPEDETAILS,
                HLATYPEDETAILS.MODIFIED,
                HLATYPEDETAILS.SAMPLEID,
                HLATYPEDETAILS.TYPE,
                HLATYPEDETAILS.REFUNIQUECOVERAGE,
                HLATYPEDETAILS.REFSHAREDCOVERAGE,
                HLATYPEDETAILS.REFWILDCARDCOVERAGE,
                HLATYPEDETAILS.TUMORCOPYNUMBER,
                HLATYPEDETAILS.TUMORUNIQUECOVERAGE,
                HLATYPEDETAILS.TUMORSHAREDCOVERAGE,
                HLATYPEDETAILS.TUMORWILDCARDCOVERAGE,
                HLATYPEDETAILS.SOMATICMISSENSE,
                HLATYPEDETAILS.SOMATICNONSENSEORFRAMESHIFT,
                HLATYPEDETAILS.SOMATICSPLICE,
                HLATYPEDETAILS.SOMATICSYNONYMOUS,
                HLATYPEDETAILS.SOMATICINFRAMEINDEL                );
        for (HlaTypeDetails type : typeDetails) {
            inserter.values(timestamp,
                    sample,
                    type.type(),
                    type.referenceUniqueCoverage(),
                    type.referenceSharedCoverage(),
                    type.referenceWildcardCoverage(),
                    type.tumorCopyNumber(),
                    type.tumorUniqueCoverage(),
                    type.tumorSharedCoverage(),
                    type.tumorWildcardCoverage(),
                    type.somaticMissense(),
                    type.somaticNonsenseOrFrameshift(),
                    type.somaticSplice(),
                    type.somaticSynonymous(),
                    type.somaticInframeIndel()
                    );
        }
        inserter.execute();
    }

    void deleteHlaFprSample(@NotNull String sample) {
        context.delete(HLATYPE).where(HLATYPE.SAMPLEID.eq(sample)).execute();
        context.delete(HLATYPEDETAILS).where(HLATYPEDETAILS.SAMPLEID.eq(sample)).execute();

    }
}