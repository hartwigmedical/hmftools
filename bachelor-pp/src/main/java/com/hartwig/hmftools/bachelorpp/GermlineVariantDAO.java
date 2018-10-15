package com.hartwig.hmftools.bachelorpp;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GERMLINEVARIANT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.bachelorpp.types.BachelorGermlineVariant;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.patientdb.dao.DatabaseUtil;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStepN;

public class GermlineVariantDAO {

    @NotNull
    private final DSLContext context;

    private static final String PASS = "PASS";

    GermlineVariantDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void write(final String sampleId, final List<BachelorGermlineVariant> bachRecords) {
        // CHSH: first remove any existing records for this patient
        context.delete(GERMLINEVARIANT).where(GERMLINEVARIANT.SAMPLEID.eq(sampleId)).execute();

        final Timestamp timestamp = new Timestamp(new Date().getTime());

        InsertValuesStepN inserter = context.insertInto(GERMLINEVARIANT,
                GERMLINEVARIANT.SAMPLEID,
                GERMLINEVARIANT.CHROMOSOME,
                GERMLINEVARIANT.POSITION,
                GERMLINEVARIANT.FILTER,
                GERMLINEVARIANT.TYPE,
                GERMLINEVARIANT.REF,
                GERMLINEVARIANT.ALT,
                GERMLINEVARIANT.GENE,
                GERMLINEVARIANT.DBSNPID,
                GERMLINEVARIANT.COSMICID,
                GERMLINEVARIANT.EFFECT,
                GERMLINEVARIANT.CODINGEFFECT,
                GERMLINEVARIANT.TRANSCRIPT,
                GERMLINEVARIANT.ALLELEREADCOUNT,
                GERMLINEVARIANT.TOTALREADCOUNT,
                GERMLINEVARIANT.ADJUSTEDCOPYNUMBER,
                GERMLINEVARIANT.ADJUSTEDVAF,
                GERMLINEVARIANT.HIGHCONFIDENCE,
                GERMLINEVARIANT.TRINUCLEOTIDECONTEXT,
                GERMLINEVARIANT.MICROHOMOLOGY,
                GERMLINEVARIANT.REPEATSEQUENCE,
                GERMLINEVARIANT.REPEATCOUNT,
                GERMLINEVARIANT.HGVSPROTEIN,
                GERMLINEVARIANT.HGVSCODING,
                GERMLINEVARIANT.BIALLELIC,
                GERMLINEVARIANT.HOTSPOT,
                GERMLINEVARIANT.MAPPABILITY,
                GERMLINEVARIANT.GERMLINESTATUS,
                GERMLINEVARIANT.MINORALLELEPLOIDY,
                GERMLINEVARIANT.REFSTATUS,
                GERMLINEVARIANT.PROGRAM,
                GERMLINEVARIANT.SOURCE,
                GERMLINEVARIANT.MODIFIED);

        for (final BachelorGermlineVariant bachRecord : bachRecords) {
            if (!bachRecord.isValid()) {
                continue;
            }

            final EnrichedSomaticVariant region = bachRecord.getEnrichedVariant();
            inserter.values(sampleId,
                    bachRecord.chromosome(),
                    bachRecord.position(),
                    bachRecord.isLowScore() ? "ARTEFACT" : PASS,
                    region.type(),
                    region.ref(),
                    region.alt(),
                    bachRecord.gene(),
                    region.dbsnpID() == null ? "" : region.dbsnpID(),
                    region.canonicalCosmicID() == null ? "" : region.canonicalCosmicID(),
                    bachRecord.effects(),
                    region.worstCodingEffect(),
                    bachRecord.transcriptId(),
                    bachRecord.getAltCount(),
                    bachRecord.getAltCount() + bachRecord.getRefCount(),
                    DatabaseUtil.decimal(region.adjustedCopyNumber()),
                    DatabaseUtil.decimal(bachRecord.getAdjustedVaf()), // region.adjustedVAF()
                    region.highConfidenceRegion(),
                    region.trinucleotideContext(),
                    region.microhomology(),
                    region.repeatSequence(),
                    region.repeatCount(),
                    bachRecord.hgvsProtein(),
                    bachRecord.hgvsCoding(),
                    bachRecord.isBiallelic(),
                    region.isHotspot(),
                    DatabaseUtil.decimal(region.mappability()),
                    region.germlineStatus(),
                    DatabaseUtil.decimal(region.minorAllelePloidy()),
                    bachRecord.isHomozygous() ? "HOM" : "HET",
                    bachRecord.program(),
                    bachRecord.source(),
                    timestamp);
        }

        inserter.execute();
    }
}
