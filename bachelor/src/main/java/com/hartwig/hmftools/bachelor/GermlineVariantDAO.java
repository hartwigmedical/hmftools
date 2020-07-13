package com.hartwig.hmftools.bachelor;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.bachelor.types.GermlineVariant;
import com.hartwig.hmftools.patientdb.dao.DatabaseUtil;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStepN;

public class GermlineVariantDAO {

    @NotNull
    private final DSLContext context;

    GermlineVariantDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void write(final String sampleId, final List<GermlineVariant> germlineVariants)
    {
        context.delete(Tables.GERMLINEVARIANT).where(Tables.GERMLINEVARIANT.SAMPLEID.eq(sampleId)).execute();

        final Timestamp timestamp = new Timestamp(new Date().getTime());

        InsertValuesStepN inserter = context.insertInto(Tables.GERMLINEVARIANT,
                Tables.GERMLINEVARIANT.SAMPLEID,
                Tables.GERMLINEVARIANT.CHROMOSOME,
                Tables.GERMLINEVARIANT.POSITION,
                Tables.GERMLINEVARIANT.FILTER,
                Tables.GERMLINEVARIANT.REPORTED,
                Tables.GERMLINEVARIANT.PATHOGENIC,
                Tables.GERMLINEVARIANT.TYPE,
                Tables.GERMLINEVARIANT.REF,
                Tables.GERMLINEVARIANT.ALT,
                Tables.GERMLINEVARIANT.GENE,
                Tables.GERMLINEVARIANT.EFFECT,
                Tables.GERMLINEVARIANT.CODINGEFFECT,
                Tables.GERMLINEVARIANT.TRANSCRIPT,
                Tables.GERMLINEVARIANT.ALLELEREADCOUNT,
                Tables.GERMLINEVARIANT.TOTALREADCOUNT,
                Tables.GERMLINEVARIANT.ADJUSTEDCOPYNUMBER,
                Tables.GERMLINEVARIANT.ADJUSTEDVAF,
                Tables.GERMLINEVARIANT.TRINUCLEOTIDECONTEXT,
                Tables.GERMLINEVARIANT.MICROHOMOLOGY,
                Tables.GERMLINEVARIANT.REPEATSEQUENCE,
                Tables.GERMLINEVARIANT.REPEATCOUNT,
                Tables.GERMLINEVARIANT.HGVSPROTEIN,
                Tables.GERMLINEVARIANT.HGVSCODING,
                Tables.GERMLINEVARIANT.BIALLELIC,
                Tables.GERMLINEVARIANT.MINORALLELECOPYNUMBER,
                Tables.GERMLINEVARIANT.REFSTATUS,
                Tables.GERMLINEVARIANT.CLINVARINFO,
                Tables.GERMLINEVARIANT.MODIFIED);

        for (final GermlineVariant variant : germlineVariants)
        {
            inserter.values(sampleId,
                    variant.chromosome(),
                    variant.position(),
                    variant.filter(),
                    variant.reported(),
                    variant.pathogenic(),
                    variant.type(),
                    variant.ref(),
                    variant.alts(),
                    variant.gene(),
                    variant.effects(),
                    variant.codingEffect(),
                    variant.transcriptId(),
                    variant.alleleReadCount(),
                    variant.totalReadCount(),
                    DatabaseUtil.decimal(variant.adjustedCopyNumber()),
                    DatabaseUtil.decimal(variant.adjustedVaf()),
                    variant.trinucleotideContext(),
                    variant.microhomology(),
                    variant.repeatSequence(),
                    variant.repeatCount(),
                    variant.hgvsProtein(),
                    variant.hgvsCoding(),
                    variant.biallelic(),
                    DatabaseUtil.decimal(variant.minorAlleleJcn()),
                    variant.refStatus(),
                    variant.clinvarInfo(),
                    timestamp);
        }

        inserter.execute();
    }
}
