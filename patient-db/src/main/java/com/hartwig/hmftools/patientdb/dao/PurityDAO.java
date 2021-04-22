package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.common.purple.purity.FittedPurityMethod.NO_TUMOR;
import static com.hartwig.hmftools.common.purple.qc.PurpleQCStatus.PASS;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PURITY;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PURITYRANGE;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.purity.ImmutableSamplePurity;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.SamplePurity;
import com.hartwig.hmftools.common.purple.qc.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep8;
import org.jooq.Record;
import org.jooq.Result;

class PurityDAO {

    @NotNull
    private final DSLContext context;

    PurityDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    @NotNull
    List<SamplePurity> readPassingQC(double minPurity) {
        List<SamplePurity> sampleIds = Lists.newArrayList();

        Result<Record> result = context.select()
                .from(PURITY)
                .where(PURITY.PURITY_.ge(minPurity))
                .and(PURITY.FITMETHOD.ne(NO_TUMOR.toString()))
                .and(PURITY.QCSTATUS.eq(PASS.toString()))
                .fetch();

        for (Record record : result) {
            sampleIds.add(ImmutableSamplePurity.builder()
                    .sampleId(record.getValue(PURITY.SAMPLEID))
                    .purity(record.getValue(PURITY.PURITY_))
                    .build());
        }

        return sampleIds;
    }

    @NotNull
    List<SamplePurity> readPassingQC(double minPurity, @NotNull String primaryTumorLocation) {
        List<SamplePurity> sampleIds = Lists.newArrayList();

        Result<Record> result = context.select()
                .from(PURITY)
                .join("clinical").on("purity.sampleId = clinical.sampleId")
                .where(PURITY.PURITY_.ge(minPurity))
                .and(PURITY.FITMETHOD.ne(NO_TUMOR.toString()))
                .and(PURITY.QCSTATUS.eq(PASS.toString()))
                .and("clinical.primaryTumorLocation = '" + primaryTumorLocation + "'")
                .fetch();

        for (Record record : result) {
            sampleIds.add(ImmutableSamplePurity.builder()
                    .sampleId(record.getValue(PURITY.SAMPLEID))
                    .purity(record.getValue(PURITY.PURITY_))
                    .build());
        }

        return sampleIds;
    }

    @Nullable
    PurityContext readPurityContext(@NotNull String sample) {
        Record result = context.select().from(PURITY).where(PURITY.SAMPLEID.eq(sample)).fetchOne();
        if (result == null) {
            return null;
        }

        final double actualPurity = result.getValue(PURITY.PURITY_);

        FittedPurity purity = ImmutableFittedPurity.builder()
                .purity(actualPurity)
                .normFactor(result.getValue(PURITY.NORMFACTOR))
                .score(result.getValue(PURITY.SCORE))
                .diploidProportion(result.getValue(PURITY.DIPLOIDPROPORTION))
                .ploidy(result.getValue(PURITY.PLOIDY))
                .somaticPenalty(result.getValue(PURITY.SOMATICPENALTY))
                .build();

        FittedPurityScore score = ImmutableFittedPurityScore.builder()
                .minPurity(result.getValue(PURITY.MINPURITY))
                .maxPurity(result.getValue(PURITY.MAXPURITY))
                .minPloidy(result.getValue(PURITY.MINPLOIDY))
                .maxPloidy(result.getValue(PURITY.MAXPLOIDY))
                .minDiploidProportion(result.getValue(PURITY.MINDIPLOIDPROPORTION))
                .maxDiploidProportion(result.getValue(PURITY.MAXDIPLOIDPROPORTION))
                .build();

        final FittedPurityMethod method = FittedPurityMethod.valueOf(result.getValue(PURITY.FITMETHOD));
        final Gender cobaltGender = Gender.valueOf(result.getValue(PURITY.GENDER));
        final Gender amberGender = Gender.valueOf(result.getValue(PURITY.AMBERGENDER));
        PurpleQC purpleQC = ImmutablePurpleQC.builder()
                .copyNumberSegments(result.getValue(PURITY.COPYNUMBERSEGMENTS))
                .unsupportedCopyNumberSegments(result.getValue(PURITY.UNSUPPORTEDCOPYNUMBERSEGMENTS))
                .purity(actualPurity)
                .amberGender(amberGender)
                .cobaltGender(cobaltGender)
                .deletedGenes(result.getValue(PURITY.DELETEDGENES))
                .contamination(result.getValue(PURITY.CONTAMINATION))
                .method(method)
                .addAllGermlineAberrations(GermlineAberration.fromString(result.getValue(PURITY.GERMLINEABERRATION)))
                .amberMeanDepth(0) //TODO: Add to database
                .build();

        return ImmutablePurityContext.builder()
                .bestFit(purity)
                .score(score)
                .qc(purpleQC)
                .wholeGenomeDuplication(result.getValue(PURITY.WHOLEGENOMEDUPLICATION) == 1)
                .microsatelliteStatus(MicrosatelliteStatus.valueOf(result.getValue(PURITY.MSSTATUS)))
                .microsatelliteIndelsPerMb(result.getValue(PURITY.MSINDELSPERMB))
                .tumorMutationalBurdenPerMb(result.getValue(PURITY.TMBPERMB))
                .tumorMutationalBurdenStatus(TumorMutationalStatus.valueOf(result.getValue(PURITY.TMBSTATUS)))
                .tumorMutationalLoad(result.getValue(PURITY.TML))
                .tumorMutationalLoadStatus(TumorMutationalStatus.valueOf(result.getValue(PURITY.TMLSTATUS)))
                .version(result.getValue(PURITY.VERSION))
                .gender(cobaltGender)
                .polyClonalProportion(result.getValue(PURITY.POLYCLONALPROPORTION))
                .method(method)
                .svTumorMutationalBurden(result.getValue(PURITY.SVTMB))
                .build();
    }

    @NotNull
    List<String> getSamplesPassingQC(double minPurity) {
        List<String> sampleIds = Lists.newArrayList();

        Result<Record> result = context.select()
                .from(PURITY)
                .where(PURITY.PURITY_.ge(minPurity))
                .and(PURITY.FITMETHOD.ne(NO_TUMOR.toString()))
                .and(PURITY.QCSTATUS.eq(PASS.toString()))
                .fetch();

        for (Record record : result) {
            sampleIds.add(record.getValue(PURITY.SAMPLEID));
        }

        return sampleIds;
    }

    @NotNull
    List<String> getSampleIds() {
        List<String> sampleIds = Lists.newArrayList();

        Result<Record> result = context.select().from(PURITY).fetch();

        for (Record record : result) {
            sampleIds.add(record.getValue(PURITY.SAMPLEID));
        }

        return sampleIds;
    }

    void write(@NotNull String sample, @NotNull PurityContext purity, @NotNull PurpleQC checks) {
        FittedPurity bestFit = purity.bestFit();
        FittedPurityScore score = purity.score();

        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(PURITY).where(PURITY.SAMPLEID.eq(sample)).execute();

        context.insertInto(PURITY,
                PURITY.VERSION,
                PURITY.SAMPLEID,
                PURITY.PURITY_,
                PURITY.GENDER,
                PURITY.FITMETHOD,
                PURITY.QCSTATUS,
                PURITY.NORMFACTOR,
                PURITY.SCORE,
                PURITY.SOMATICPENALTY,
                PURITY.PLOIDY,
                PURITY.DIPLOIDPROPORTION,
                PURITY.MINDIPLOIDPROPORTION,
                PURITY.MAXDIPLOIDPROPORTION,
                PURITY.MINPURITY,
                PURITY.MAXPURITY,
                PURITY.MINPLOIDY,
                PURITY.MAXPLOIDY,
                PURITY.POLYCLONALPROPORTION,
                PURITY.WHOLEGENOMEDUPLICATION,
                PURITY.MSINDELSPERMB,
                PURITY.MSSTATUS,
                PURITY.TMBPERMB,
                PURITY.TMBSTATUS,
                PURITY.TML,
                PURITY.TMLSTATUS,
                PURITY.SVTMB,
                PURITY.DELETEDGENES,
                PURITY.COPYNUMBERSEGMENTS,
                PURITY.UNSUPPORTEDCOPYNUMBERSEGMENTS,
                PURITY.CONTAMINATION,
                PURITY.GERMLINEABERRATION,
                PURITY.AMBERGENDER,
                PURITY.MODIFIED)
                .values(purity.version(),
                        sample,
                        DatabaseUtil.decimal(bestFit.purity()),
                        purity.gender().toString(),
                        purity.method().toString(),
                        checks.toString(),
                        DatabaseUtil.decimal(bestFit.normFactor()),
                        DatabaseUtil.decimal(bestFit.score()),
                        DatabaseUtil.decimal(bestFit.somaticPenalty()),
                        DatabaseUtil.decimal(bestFit.ploidy()),
                        DatabaseUtil.decimal(bestFit.diploidProportion()),
                        DatabaseUtil.decimal(score.minDiploidProportion()),
                        DatabaseUtil.decimal(score.maxDiploidProportion()),
                        DatabaseUtil.decimal(score.minPurity()),
                        DatabaseUtil.decimal(score.maxPurity()),
                        DatabaseUtil.decimal(score.minPloidy()),
                        DatabaseUtil.decimal(score.maxPloidy()),
                        DatabaseUtil.decimal(purity.polyClonalProportion()),
                        purity.wholeGenomeDuplication() ? (byte) 1 : (byte) 0,
                        DatabaseUtil.decimal(purity.microsatelliteIndelsPerMb()),
                        purity.microsatelliteStatus().toString(),
                        DatabaseUtil.decimal(purity.tumorMutationalBurdenPerMb()),
                        purity.tumorMutationalBurdenStatus().toString(),
                        DatabaseUtil.decimal(purity.tumorMutationalLoad()),
                        purity.tumorMutationalLoadStatus().toString(),
                        purity.svTumorMutationalBurden(),
                        purity.qc().deletedGenes(),
                        purity.qc().copyNumberSegments(),
                        purity.qc().unsupportedCopyNumberSegments(),
                        purity.qc().contamination(),
                        GermlineAberration.toString(purity.qc().germlineAberrations()),
                        purity.qc().amberGender(),
                        timestamp)
                .execute();
    }

    void write(@NotNull String sample, @NotNull List<FittedPurity> purities) {
        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(PURITYRANGE).where(PURITYRANGE.SAMPLEID.eq(sample)).execute();

        InsertValuesStep8 inserter = context.insertInto(PURITYRANGE,
                PURITYRANGE.SAMPLEID,
                PURITYRANGE.PURITY,
                PURITYRANGE.NORMFACTOR,
                PURITYRANGE.SCORE,
                PURITYRANGE.SOMATICPENALTY,
                PURITYRANGE.PLOIDY,
                PURITYRANGE.DIPLOIDPROPORTION,
                PURITYRANGE.MODIFIED);

        purities.forEach(x -> addPurity(timestamp, inserter, sample, x));
        inserter.execute();
    }

    private static void addPurity(@NotNull Timestamp timestamp, @NotNull InsertValuesStep8 inserter, @NotNull String sample,
            @NotNull FittedPurity purity) {
        inserter.values(sample,
                DatabaseUtil.decimal(purity.purity()),
                DatabaseUtil.decimal(purity.normFactor()),
                DatabaseUtil.decimal(purity.score()),
                DatabaseUtil.decimal(purity.somaticPenalty()),
                DatabaseUtil.decimal(purity.ploidy()),
                DatabaseUtil.decimal(purity.diploidProportion()),
                timestamp);
    }

    void deletePurityForSample(@NotNull String sample) {
        context.delete(PURITY).where(PURITY.SAMPLEID.eq(sample)).execute();
        context.delete(PURITYRANGE).where(PURITYRANGE.SAMPLEID.eq(sample)).execute();
    }
}
