package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.common.purple.FittedPurityMethod.NO_TUMOR;
import static com.hartwig.hmftools.common.purple.PurpleQCStatus.PASS;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PURITY;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PURITYRANGE;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.FittedPurityScore;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.RunMode;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep8;
import org.jooq.Record;
import org.jooq.Result;

class PurityDAO
{

    private final DSLContext context;

    PurityDAO(@NotNull final DSLContext context)
    {
        this.context = context;
    }

    @Nullable
    PurityContext readPurityContext(final String sample)
    {
        Record result = context.select().from(PURITY).where(PURITY.SAMPLEID.eq(sample)).fetchOne();
        if(result == null)
        {
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

        String runModeStr = result.get(PURITY.RUNMODE);

        RunMode runMode = runModeStr != null && !runModeStr.isEmpty() ? RunMode.valueOf(runModeStr) : RunMode.TUMOR_GERMLINE;

        Set<PurpleQCStatus> status = PurpleQCStatus.fromString(result.getValue(PURITY.QCSTATUS));

        PurpleQC purpleQC = ImmutablePurpleQC.builder()
                .status(status)
                .copyNumberSegments(result.getValue(PURITY.COPYNUMBERSEGMENTS))
                .unsupportedCopyNumberSegments(result.getValue(PURITY.UNSUPPORTEDCOPYNUMBERSEGMENTS))
                .purity(actualPurity)
                .amberGender(amberGender)
                .cobaltGender(cobaltGender)
                .deletedGenes(result.getValue(PURITY.DELETEDGENES))
                .contamination(result.getValue(PURITY.CONTAMINATION))
                .method(method)
                .addAllGermlineAberrations(GermlineAberration.fromString(result.getValue(PURITY.GERMLINEABERRATION)))
                .amberMeanDepth(0) // not persisted to database
                .lohPercent(0)
                .build();

        return ImmutablePurityContext.builder()
                .bestFit(purity)
                .score(score)
                .qc(purpleQC)
                .runMode(runMode)
                .wholeGenomeDuplication(result.getValue(PURITY.WHOLEGENOMEDUPLICATION) == 1)
                .microsatelliteStatus(MicrosatelliteStatus.valueOf(result.getValue(PURITY.MSSTATUS)))
                .microsatelliteIndelsPerMb(result.getValue(PURITY.MSINDELSPERMB))
                .tumorMutationalBurdenPerMb(result.getValue(PURITY.TMBPERMB))
                .tumorMutationalBurdenStatus(TumorMutationalStatus.valueOf(result.getValue(PURITY.TMBSTATUS)))
                .tumorMutationalLoad(result.getValue(PURITY.TML))
                .tumorMutationalLoadStatus(TumorMutationalStatus.valueOf(result.getValue(PURITY.TMLSTATUS)))
                .gender(cobaltGender)
                .polyClonalProportion(result.getValue(PURITY.POLYCLONALPROPORTION))
                .method(method)
                .svTumorMutationalBurden(result.getValue(PURITY.SVTMB))
                .targeted(result.getValue(PURITY.TARGETED) == 1)
                .build();
    }

    List<String> getSamplesPassingQC(double minPurity)
    {
        List<String> sampleIds = Lists.newArrayList();

        Result<Record> result = context.select()
                .from(PURITY)
                .where(PURITY.PURITY_.ge(minPurity))
                .and(PURITY.FITMETHOD.ne(NO_TUMOR.toString()))
                .and(PURITY.QCSTATUS.eq(PASS.toString()))
                .fetch();

        for(Record record : result)
        {
            sampleIds.add(record.getValue(PURITY.SAMPLEID));
        }

        return sampleIds;
    }

    List<String> getSampleIds()
    {
        List<String> sampleIds = Lists.newArrayList();

        Result<Record> result = context.select().from(PURITY).fetch();

        for(Record record : result)
        {
            sampleIds.add(record.getValue(PURITY.SAMPLEID));
        }

        return sampleIds;
    }

    void write(final String sample, final PurityContext purity, final PurpleQC checks)
    {
        FittedPurity bestFit = purity.bestFit();
        FittedPurityScore score = purity.score();

        Timestamp timestamp = new Timestamp(new Date().getTime());
        context.delete(PURITY).where(PURITY.SAMPLEID.eq(sample)).execute();

        context.insertInto(PURITY,
                        PURITY.SAMPLEID,
                        PURITY.PURITY_,
                        PURITY.GENDER,
                        PURITY.FITMETHOD,
                        PURITY.QCSTATUS,
                        PURITY.RUNMODE,
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
                        PURITY.TARGETED,
                        PURITY.MODIFIED)
                .values(sample,
                        DatabaseUtil.decimal(bestFit.purity()),
                        purity.gender().toString(),
                        purity.method().toString(),
                        DatabaseUtil.checkStringLength(checks.toString(), PURITY.QCSTATUS),
                        purity.runMode().toString(),
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
                        purity.targeted(),
                        timestamp)
                .execute();
    }

    void write(final String sample, final List<FittedPurity> purities)
    {
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

    private static void addPurity(final Timestamp timestamp, final InsertValuesStep8 inserter, final String sample,
            final FittedPurity purity)
    {
        inserter.values(sample,
                DatabaseUtil.decimal(purity.purity()),
                DatabaseUtil.decimal(purity.normFactor()),
                DatabaseUtil.decimal(purity.score()),
                DatabaseUtil.decimal(purity.somaticPenalty()),
                DatabaseUtil.decimal(purity.ploidy()),
                DatabaseUtil.decimal(purity.diploidProportion()),
                timestamp);
    }

    void deletePurityForSample(final String sample)
    {
        context.delete(PURITY).where(PURITY.SAMPLEID.eq(sample)).execute();
        context.delete(PURITYRANGE).where(PURITYRANGE.SAMPLEID.eq(sample)).execute();
    }
}
