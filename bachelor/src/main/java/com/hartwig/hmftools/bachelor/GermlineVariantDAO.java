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
                Tables.GERMLINEVARIANT.REFSTATUS,
                Tables.GERMLINEVARIANT.REPORTED,
                Tables.GERMLINEVARIANT.PATHOGENIC,
                Tables.GERMLINEVARIANT.CLINVARINFO,
                Tables.GERMLINEVARIANT.TYPE,
                Tables.GERMLINEVARIANT.REF,
                Tables.GERMLINEVARIANT.ALT,
                Tables.GERMLINEVARIANT.GENE,
                Tables.GERMLINEVARIANT.TRANSCRIPT,
                Tables.GERMLINEVARIANT.EFFECT,
                Tables.GERMLINEVARIANT.CODINGEFFECT,
                Tables.GERMLINEVARIANT.HGVSCODING,
                Tables.GERMLINEVARIANT.HGVSPROTEIN,
                Tables.GERMLINEVARIANT.MICROHOMOLOGY,
                Tables.GERMLINEVARIANT.REPEATSEQUENCE,
                Tables.GERMLINEVARIANT.REPEATCOUNT,
                Tables.GERMLINEVARIANT.TRINUCLEOTIDECONTEXT,
                Tables.GERMLINEVARIANT.ALLELEREADCOUNT,
                Tables.GERMLINEVARIANT.TOTALREADCOUNT,
                Tables.GERMLINEVARIANT.ADJUSTEDCOPYNUMBER,
                Tables.GERMLINEVARIANT.MINORALLELECOPYNUMBER,
                Tables.GERMLINEVARIANT.ADJUSTEDVAF,
                Tables.GERMLINEVARIANT.BIALLELIC,
                Tables.GERMLINEVARIANT.MODIFIED);

        for (final GermlineVariant variant : germlineVariants)
        {
            inserter.values(sampleId,
                    variant.chromosome(),
                    variant.position(),
                    variant.filter(),
                    variant.refStatus(),
                    variant.reported(),
                    variant.pathogenic(),
                    variant.clinvarInfo(),
                    variant.type(),
                    variant.ref(),
                    variant.alts(),
                    variant.gene(),
                    variant.transcriptId(),
                    variant.effects(),
                    variant.codingEffect(),
                    variant.hgvsCoding(),
                    variant.hgvsProtein(),
                    variant.microhomology(),
                    variant.repeatSequence(),
                    variant.repeatCount(),
                    variant.trinucleotideContext(),
                    variant.alleleReadCount(),
                    variant.totalReadCount(),
                    DatabaseUtil.decimal(variant.adjustedCopyNumber()),
                    DatabaseUtil.decimal(variant.minorAlleleJcn()),
                    DatabaseUtil.decimal(variant.adjustedVaf()),
                    variant.biallelic(),
                    timestamp);
        }

        inserter.execute();
    }
}
