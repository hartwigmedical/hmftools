package com.hartwig.hmftools.common.purple;

import java.io.IOException;

public final class PurityContextFile
{
    public static PurityContext read(final String basePath, final String sample) throws IOException
    {
        return readWithQC(PurpleQCFile.generateFilename(basePath, sample), PurplePurity.generateFilename(basePath, sample));
    }

    public static void write(final String basePath, final String sample, final PurityContext context) throws IOException
    {
        PurpleQCFile.write(PurpleQCFile.generateFilename(basePath, sample), context.qc());
        PurplePurity.write(PurplePurity.generateFilename(basePath, sample), context);
    }

    public static PurityContext readWithQC(final String qcFile, final String purityFile) throws IOException
    {
        PurpleQC qc = PurpleQCFile.read(qcFile);

        PurplePurity purplePurity = PurplePurity.read(purityFile);

        ImmutablePurityContext.Builder builder = ImmutablePurityContext.builder();

        builder.qc(qc);

        FittedPurity fittedPurity = ImmutableFittedPurity.builder()
                .purity(purplePurity.Purity)
                .normFactor(purplePurity.NormFactor)
                .score(purplePurity.Score)
                .diploidProportion(purplePurity.DiploidProportion)
                .ploidy(purplePurity.Ploidy)
                .somaticPenalty(purplePurity.SomaticPenalty)
                .build();

        FittedPurityScore score = ImmutableFittedPurityScore.builder()
                .minPurity(purplePurity.MinPurity)
                .maxPurity(purplePurity.MaxPurity)
                .minPloidy(purplePurity.MinPloidy)
                .maxPloidy(purplePurity.MaxPloidy)
                .minDiploidProportion(purplePurity.MinDiploidProportion)
                .maxDiploidProportion(purplePurity.MaxDiploidProportion)
                .build();

        builder.score(score)
                .bestFit(fittedPurity)
                .gender(purplePurity.Sex)
                .method(purplePurity.FitMethod)
                .runMode(purplePurity.Mode)
                .targeted(purplePurity.Targeted)
                .polyClonalProportion(purplePurity.PolyclonalProportion)
                .wholeGenomeDuplication(purplePurity.WholeGenomeDuplication)
                .microsatelliteIndelsPerMb(purplePurity.MsIndelsPerMb)
                .microsatelliteStatus(purplePurity.MsStatus)
                .tumorMutationalLoad(purplePurity.Tml)
                .tumorMutationalLoadStatus(purplePurity.TmlStatus)
                .tumorMutationalBurdenPerMb(purplePurity.TmbPerMb)
                .tumorMutationalBurdenStatus(purplePurity.TmbStatus)
                .svTumorMutationalBurden(purplePurity.SvTumorMutationalBurden);

        return builder.build();
    }
}
