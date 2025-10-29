package com.hartwig.hmftools.compar.purple;

import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.FittedPurityScore;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.MicrosatelliteStatus;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.RunMode;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestPurityDataBuilder
{
    public double purity = 0.8;
    public double ploidy = 2.1;
    public double contamination = 0;
    public double tmb = 2;
    public double msIndels = 0.1;
    public int tml = 50;
    public int copyNumberSegments = 60;
    public int unsupportedCopyNumberSegments = 2;
    public int svTmb = 5;
    public Set<PurpleQCStatus> qcStatus = Set.of(PurpleQCStatus.PASS);
    public Gender gender = Gender.FEMALE;
    public Set<GermlineAberration> germlineAberrations = Set.of(GermlineAberration.NONE);
    public FittedPurityMethod fitMethod = FittedPurityMethod.SOMATIC;
    public MicrosatelliteStatus msStatus = MicrosatelliteStatus.MSS;
    public TumorMutationalStatus tmbStatus = TumorMutationalStatus.LOW;
    public TumorMutationalStatus tmlStatus = TumorMutationalStatus.LOW;
    public double tincLevel = 0;

    private static final Consumer<TestPurityDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.purity = 0.2;
        b.ploidy = 3.1;
        b.contamination = 0.1;
        b.tmb = 20;
        b.msIndels = 10;
        b.tml = 500;
        b.copyNumberSegments = 600;
        b.unsupportedCopyNumberSegments = 20;
        b.svTmb = 50;
        b.qcStatus = Set.of(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE, PurpleQCStatus.WARN_DELETED_GENES);
        b.gender = Gender.MALE;
        b.germlineAberrations = Set.of(GermlineAberration.KLINEFELTER);
        b.fitMethod = FittedPurityMethod.NORMAL;
        b.msStatus = MicrosatelliteStatus.MSI;
        b.tmbStatus = TumorMutationalStatus.HIGH;
        b.tmlStatus = TumorMutationalStatus.HIGH;
        b.tincLevel = 0.3;
    };

    public static final TestComparableItemBuilder<TestPurityDataBuilder, PurityData> BUILDER =
            new TestComparableItemBuilder<>(TestPurityDataBuilder::new, TestPurityDataBuilder::build, ALTERNATE_INITIALIZER);

    private PurityData build()
    {
        FittedPurity bestFit = ImmutableFittedPurity.builder()
                .purity(purity)
                .ploidy(ploidy)
                .normFactor(-1)
                .score(-1)
                .diploidProportion(-1)
                .somaticPenalty(-1)
                .build();
        PurpleQC qc = ImmutablePurpleQC.builder()
                .status(qcStatus)
                .method(fitMethod)
                .copyNumberSegments(copyNumberSegments)
                .unsupportedCopyNumberSegments(unsupportedCopyNumberSegments)
                .purity(purity)
                .contamination(contamination)
                .cobaltGender(gender)
                .amberGender(gender)
                .germlineAberrations(germlineAberrations)
                .tincLevel(tincLevel)
                .deletedGenes(-1)
                .amberMeanDepth(-1)
                .lohPercent(-1)
                .chimerismPercentage(-1)
                .build();
        FittedPurityScore score = ImmutableFittedPurityScore.builder()
                .minPurity(-1)
                .maxPurity(-1)
                .minPloidy(-1)
                .maxPloidy(-1)
                .minDiploidProportion(-1)
                .maxDiploidProportion(-1)
                .build();
        return new PurityData(
                ImmutablePurityContext.builder()
                        .gender(gender)
                        .bestFit(bestFit)
                        .method(fitMethod)
                        .qc(qc)
                        .microsatelliteIndelsPerMb(msIndels)
                        .tumorMutationalBurdenPerMb(tmb)
                        .tumorMutationalLoad(tml)
                        .svTumorMutationalBurden(svTmb)
                        .microsatelliteStatus(msStatus)
                        .tumorMutationalLoadStatus(tmlStatus)
                        .tumorMutationalBurdenStatus(tmbStatus)
                        .runMode(RunMode.TUMOR_GERMLINE)
                        .targeted(false)
                        .score(score)
                        .polyClonalProportion(-1)
                        .wholeGenomeDuplication(false)
                        .build()
        );
    }
}
