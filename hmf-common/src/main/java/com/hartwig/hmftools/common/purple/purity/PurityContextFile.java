package com.hartwig.hmftools.common.purple.purity;

import static com.hartwig.hmftools.common.purple.PurpleCommon.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;

import org.jetbrains.annotations.NotNull;

public final class PurityContextFile
{
    private static final DecimalFormat FORMAT = PurpleCommon.decimalFormat("0.0000");

    private static final String EXTENSION = ".purple.purity.tsv";

    public static PurityContext read(final String basePath, final String sample) throws IOException
    {
        return readWithQC(PurpleQCFile.generateFilename(basePath, sample), generateFilenameForReading(basePath, sample));
    }

    public static String generateFilenameForReading(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(final String basePath, final String sample, final PurityContext context) throws IOException
    {
        PurpleQCFile.write(PurpleQCFile.generateFilename(basePath, sample), context.qc());
        writeBestPurity(basePath, sample, context);
    }

    private static void writeBestPurity(final String basePath, final String sample, final PurityContext context) throws IOException
    {
        final String filePath = generateFilenameForWriting(basePath, sample);
        Files.write(new File(filePath).toPath(), toLines(context));
    }

    @NotNull
    private static String generateFilenameForWriting(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    @VisibleForTesting
    static List<String> toLines(final PurityContext context)
    {
        return Lists.newArrayList(header(), toString(context));
    }

    @NotNull
    static String header()
    {
        return new StringJoiner(DELIMITER, "", "")
                .add("purity")
                .add("normFactor")
                .add("score")
                .add("diploidProportion")
                .add("ploidy")
                .add("gender")
                .add("status")
                .add("runMode")
                .add("targeted")
                .add("polyclonalProportion")
                .add("minPurity")
                .add("maxPurity")
                .add("minPloidy")
                .add("maxPloidy")
                .add("minDiploidProportion")
                .add("maxDiploidProportion")
                .add("version")
                .add("somaticPenalty")
                .add("wholeGenomeDuplication")
                .add("msIndelsPerMb")
                .add("msStatus")
                .add("tml")
                .add("tmlStatus")
                .add("tmbPerMb")
                .add("tmbStatus")
                .add("svTumorMutationalBurden")
                .toString();
    }

    @NotNull
    static String toString(final PurityContext context)
    {
        final FittedPurity purity = context.bestFit();
        final FittedPurityScore score = context.score();

        return new StringJoiner(DELIMITER)
                .add(FORMAT.format(purity.purity()))
                .add(FORMAT.format(purity.normFactor()))
                .add(FORMAT.format(purity.score()))
                .add(FORMAT.format(purity.diploidProportion()))
                .add(FORMAT.format(purity.ploidy()))
                .add(String.valueOf(context.gender()))
                .add(String.valueOf(context.method()))
                .add(String.valueOf(context.runMode()))
                .add(String.valueOf(context.targeted()))
                .add(FORMAT.format(context.polyClonalProportion()))
                .add(FORMAT.format(score.minPurity()))
                .add(FORMAT.format(score.maxPurity()))
                .add(FORMAT.format(score.minPloidy()))
                .add(FORMAT.format(score.maxPloidy()))
                .add(FORMAT.format(score.minDiploidProportion()))
                .add(FORMAT.format(score.maxDiploidProportion()))
                .add(String.valueOf(context.version()))
                .add(FORMAT.format(purity.somaticPenalty()))
                .add(String.valueOf(context.wholeGenomeDuplication()))
                .add(FORMAT.format(context.microsatelliteIndelsPerMb()))
                .add(String.valueOf(context.microsatelliteStatus()))
                .add(String.valueOf(context.tumorMutationalLoad()))
                .add(String.valueOf(context.tumorMutationalLoadStatus()))
                .add(FORMAT.format(context.tumorMutationalBurdenPerMb()))
                .add(String.valueOf(context.tumorMutationalBurdenStatus()))
                .add(String.valueOf(context.svTumorMutationalBurden()))
                .toString();
    }

    public static PurityContext readWithQC(final String qcFilePath, String filePath) throws IOException
    {
        PurpleQC qc = PurpleQCFile.read(qcFilePath);

        List<String> lines = Files.readAllLines(new File(filePath).toPath());

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
        String[] values = lines.get(1).split(DELIMITER, -1);

        ImmutablePurityContext.Builder builder = ImmutablePurityContext.builder();

        builder.qc(qc);

        FittedPurity fittedPurity = ImmutableFittedPurity.builder()
                .purity(Double.parseDouble(values[fieldsIndexMap.get("purity")]))
                .normFactor(Double.parseDouble(values[fieldsIndexMap.get("normFactor")]))
                .score(Double.parseDouble(values[fieldsIndexMap.get("score")]))
                .diploidProportion(Double.parseDouble(values[fieldsIndexMap.get("diploidProportion")]))
                .ploidy(Double.parseDouble(values[fieldsIndexMap.get("ploidy")]))
                .somaticPenalty(Double.parseDouble(values[fieldsIndexMap.get("somaticPenalty")]))
                .build();

        FittedPurityScore score = ImmutableFittedPurityScore.builder()
                .minPurity(Double.parseDouble(values[fieldsIndexMap.get("minPurity")]))
                .maxPurity(Double.parseDouble(values[fieldsIndexMap.get("maxPurity")]))
                .minPloidy(Double.parseDouble(values[fieldsIndexMap.get("minPloidy")]))
                .maxPloidy(Double.parseDouble(values[fieldsIndexMap.get("maxPloidy")]))
                .minDiploidProportion(Double.parseDouble(values[fieldsIndexMap.get("minDiploidProportion")]))
                .maxDiploidProportion(Double.parseDouble(values[fieldsIndexMap.get("maxDiploidProportion")]))
                .build();

        String tml =  values[fieldsIndexMap.get("tml")];
        int tumorMutationalLoad = tml.contains(".") ?
                Integer.parseInt(tml.substring(0, tml.indexOf("."))) : Integer.parseInt(tml);

        RunMode runMode = fieldsIndexMap.containsKey("runMode") ?
                RunMode.valueOf(values[fieldsIndexMap.get("runMode")]) : RunMode.TUMOR_GERMLINE;

        builder.score(score)
                .bestFit(fittedPurity)
                .gender(Gender.valueOf(values[fieldsIndexMap.get("gender")]))
                .method(FittedPurityMethod.valueOf(values[fieldsIndexMap.get("status")]))
                .runMode(runMode)
                .targeted(fieldsIndexMap.containsKey("targeted") ? Boolean.parseBoolean(values[fieldsIndexMap.get("targeted")]) : false)
                .polyClonalProportion(Double.parseDouble(values[fieldsIndexMap.get("polyclonalProportion")]))
                .version(values[fieldsIndexMap.get("version")])
                .wholeGenomeDuplication(Boolean.parseBoolean(values[fieldsIndexMap.get("wholeGenomeDuplication")]))
                .microsatelliteIndelsPerMb(Double.parseDouble(values[fieldsIndexMap.get("msIndelsPerMb")]))
                .microsatelliteStatus(MicrosatelliteStatus.valueOf(values[fieldsIndexMap.get("msStatus")]))
                .tumorMutationalLoad(tumorMutationalLoad)
                .tumorMutationalLoadStatus(TumorMutationalStatus.valueOf(values[fieldsIndexMap.get("tmlStatus")]))
                .tumorMutationalBurdenPerMb(Double.parseDouble(values[fieldsIndexMap.get("tmbPerMb")]))
                .tumorMutationalBurdenStatus(TumorMutationalStatus.valueOf(values[fieldsIndexMap.get("tmbStatus")]))
                .svTumorMutationalBurden(Integer.parseInt(values[fieldsIndexMap.get("svTumorMutationalBurden")]));

        return builder.build();
    }
}
