package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public final class PurplePurity
{
    public final double Purity;
    public final double NormFactor;
    public final double Score;
    public final double DiploidProportion;
    public final double Ploidy;
    public final Gender Sex;
    public final FittedPurityMethod FitMethod;
    public final double PolyclonalProportion;
    public final double MinPurity;
    public final double MaxPurity;
    public final double MinPloidy;
    public final double MaxPloidy;
    public final double MinDiploidProportion;
    public final double MaxDiploidProportion;
    public final double SomaticPenalty;
    public final boolean WholeGenomeDuplication;
    public final double MsIndelsPerMb;
    public final MicrosatelliteStatus MsStatus;
    public final int Tml;
    public final TumorMutationalStatus TmlStatus;
    public final double TmbPerMb;
    public final TumorMutationalStatus TmbStatus;
    public final int SvTumorMutationalBurden;
    public final RunMode Mode;
    public final boolean Targeted;

    private static final DecimalFormat FORMAT = PurpleCommon.decimalFormat("0.0000");

    private static final String EXTENSION = ".purple.purity.tsv";

    public PurplePurity(
            final double purity, final double normFactor, final double score, final double diploidProportion,
            final double ploidy, final Gender sex, final FittedPurityMethod fitMethod, final double polyclonalProportion,
            final double minPurity, final double maxPurity, final double minPloidy, final double maxPloidy,
            final double minDiploidProportion, final double maxDiploidProportion, final double somaticPenalty,
            final boolean wholeGenomeDuplication, final double msIndelsPerMb, final MicrosatelliteStatus msStatus, final int tml,
            final TumorMutationalStatus tmlStatus, final double tmbPerMb, final TumorMutationalStatus tmbStatus,
            final int svTumorMutationalBurden, final RunMode mode, final boolean targeted)
    {
        Purity = purity;
        NormFactor = normFactor;
        Score = score;
        DiploidProportion = diploidProportion;
        Ploidy = ploidy;
        Sex = sex;
        FitMethod = fitMethod;
        PolyclonalProportion = polyclonalProportion;
        MinPurity = minPurity;
        MaxPurity = maxPurity;
        MinPloidy = minPloidy;
        MaxPloidy = maxPloidy;
        MinDiploidProportion = minDiploidProportion;
        MaxDiploidProportion = maxDiploidProportion;
        SomaticPenalty = somaticPenalty;
        WholeGenomeDuplication = wholeGenomeDuplication;
        MsIndelsPerMb = msIndelsPerMb;
        MsStatus = msStatus;
        Tml = tml;
        TmlStatus = tmlStatus;
        TmbPerMb = tmbPerMb;
        TmbStatus = tmbStatus;
        SvTumorMutationalBurden = svTumorMutationalBurden;
        Mode = mode;
        Targeted = targeted;
    }

    public static PurplePurity read(final String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
        String[] values = lines.get(1).split(TSV_DELIM, -1);

        return new PurplePurity(
                Double.parseDouble(values[fieldsIndexMap.get(PurplePurityColumn.purity.toString())]),
                Double.parseDouble(values[fieldsIndexMap.get(PurplePurityColumn.normFactor.toString())]),
                Double.parseDouble(values[fieldsIndexMap.get(PurplePurityColumn.score.toString())]),
                Double.parseDouble(values[fieldsIndexMap.get(PurplePurityColumn.diploidProportion.toString())]),
                Double.parseDouble(values[fieldsIndexMap.get(PurplePurityColumn.ploidy.toString())]),
                Gender.valueOf(values[fieldsIndexMap.get(PurplePurityColumn.gender.toString())]),
                FittedPurityMethod.valueOf(values[fieldsIndexMap.get(PurplePurityColumn.status.toString())]),
                Double.parseDouble(values[fieldsIndexMap.get(PurplePurityColumn.polyclonalProportion.toString())]),
                Double.parseDouble(values[fieldsIndexMap.get(PurplePurityColumn.minPurity.toString())]),
                Double.parseDouble(values[fieldsIndexMap.get(PurplePurityColumn.maxPurity.toString())]),
                Double.parseDouble(values[fieldsIndexMap.get(PurplePurityColumn.minPloidy.toString())]),
                Double.parseDouble(values[fieldsIndexMap.get(PurplePurityColumn.maxPloidy.toString())]),
                Double.parseDouble(values[fieldsIndexMap.get(PurplePurityColumn.minDiploidProportion.toString())]),
                Double.parseDouble(values[fieldsIndexMap.get(PurplePurityColumn.maxDiploidProportion.toString())]),
                Double.parseDouble(values[fieldsIndexMap.get(PurplePurityColumn.somaticPenalty.toString())]),
                Boolean.parseBoolean(values[fieldsIndexMap.get(PurplePurityColumn.wholeGenomeDuplication.toString())]),
                Double.parseDouble(values[fieldsIndexMap.get(PurplePurityColumn.msIndelsPerMb.toString())]),
                MicrosatelliteStatus.valueOf(values[fieldsIndexMap.get(PurplePurityColumn.msStatus.toString())]),
                Integer.parseInt(values[fieldsIndexMap.get(PurplePurityColumn.tml.toString())]),
                TumorMutationalStatus.valueOf(values[fieldsIndexMap.get(PurplePurityColumn.tmbStatus.toString())]),
                Double.parseDouble(values[fieldsIndexMap.get(PurplePurityColumn.tmbPerMb.toString())]),
                TumorMutationalStatus.valueOf(values[fieldsIndexMap.get(PurplePurityColumn.tmlStatus.toString())]),
                Integer.parseInt(values[fieldsIndexMap.get(PurplePurityColumn.svTumorMutationalBurden.toString())]),
                RunMode.valueOf(values[fieldsIndexMap.get(PurplePurityColumn.runMode.toString())]),
                Boolean.parseBoolean(values[fieldsIndexMap.get(PurplePurityColumn.targeted.toString())]));
    }

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(final String filename, final PurplePurity purplePurity) throws IOException
    {
        List<String> lines = List.of(header(), toString(purplePurity));
        Files.write(new File(filename).toPath(), lines);
    }

    public static void write(final String filename, final PurityContext context) throws IOException
    {
        PurplePurity purplePurity = new PurplePurity(
                context.bestFit().purity(), context.bestFit().normFactor(), context.bestFit().score(), context.bestFit().diploidProportion(),
                context.bestFit().ploidy(), context.gender(), context.method(), context.polyClonalProportion(),
                context.score().minPurity(), context.score().maxPurity(), context.score().minPloidy(), context.score().maxPloidy(),
                context.score().minDiploidProportion(), context.score().maxDiploidProportion(), context.bestFit().somaticPenalty(),
                context.wholeGenomeDuplication(), context.microsatelliteIndelsPerMb(), context.microsatelliteStatus(),
                context.tumorMutationalLoad(), context.tumorMutationalLoadStatus(), context.tumorMutationalBurdenPerMb(),
                context.tumorMutationalBurdenStatus(), context.svTumorMutationalBurden(), context.runMode(), context.targeted());

        write(filename, purplePurity);
    }

    private enum PurplePurityColumn
    {
        purity,
        normFactor,
        score,
        diploidProportion,
        ploidy,
        gender,
        status,
        polyclonalProportion,
        minPurity,
        maxPurity,
        minPloidy,
        maxPloidy,
        minDiploidProportion,
        maxDiploidProportion,
        somaticPenalty,
        wholeGenomeDuplication,
        msIndelsPerMb,
        msStatus,
        tml,
        tmlStatus,
        tmbPerMb,
        tmbStatus,
        svTumorMutationalBurden,
        runMode,
        targeted;
    }

    private static String header()
    {
        return Arrays.stream(PurplePurityColumn.values()).map(x -> x.toString()).collect(Collectors.joining(TSV_DELIM));
    }

    private static String toString(final PurplePurity purplePurity)
    {
        return new StringJoiner(TSV_DELIM)
                .add(FORMAT.format(purplePurity.Purity))
                .add(FORMAT.format(purplePurity.NormFactor))
                .add(FORMAT.format(purplePurity.Score))
                .add(FORMAT.format(purplePurity.DiploidProportion))
                .add(FORMAT.format(purplePurity.Ploidy))
                .add(String.valueOf(purplePurity.Sex))
                .add(String.valueOf(purplePurity.FitMethod))
                .add(FORMAT.format(purplePurity.PolyclonalProportion))
                .add(FORMAT.format(purplePurity.MinPurity))
                .add(FORMAT.format(purplePurity.MaxPurity))
                .add(FORMAT.format(purplePurity.MinPloidy))
                .add(FORMAT.format(purplePurity.MaxPloidy))
                .add(FORMAT.format(purplePurity.MinDiploidProportion))
                .add(FORMAT.format(purplePurity.MaxDiploidProportion))
                .add(FORMAT.format(purplePurity.SomaticPenalty))
                .add(String.valueOf(purplePurity.WholeGenomeDuplication))
                .add(FORMAT.format(purplePurity.MsIndelsPerMb))
                .add(String.valueOf(purplePurity.MsStatus))
                .add(String.valueOf(purplePurity.Tml))
                .add(String.valueOf(purplePurity.TmbStatus))
                .add(FORMAT.format(purplePurity.TmbPerMb))
                .add(String.valueOf(purplePurity.TmbStatus))
                .add(String.valueOf(purplePurity.SvTumorMutationalBurden))
                .add(String.valueOf(purplePurity.Mode))
                .add(String.valueOf(purplePurity.Targeted))
                .toString();
    }
}
