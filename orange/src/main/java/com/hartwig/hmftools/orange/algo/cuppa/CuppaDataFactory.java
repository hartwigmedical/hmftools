package com.hartwig.hmftools.orange.algo.cuppa;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.cuppa.ClassifierType;
import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
import com.hartwig.hmftools.common.cuppa.DataTypes;
import com.hartwig.hmftools.common.cuppa.ResultType;
import com.hartwig.hmftools.common.cuppa.SvDataType;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaData;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaPrediction;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CuppaDataFactory
{
    @NotNull
    public static CuppaData create(@NotNull List<CuppaDataFile> entries)
    {
        return ImmutableCuppaData.builder()
                .predictions(extractPredictions(entries))
                .simpleDups32To200B(safeInt(entries, SvDataType.SIMPLE_DUP_32B_200B))
                .maxComplexSize(safeInt(entries, SvDataType.MAX_COMPLEX_SIZE))
                .lineCount(safeInt(entries, SvDataType.LINE))
                .telomericSGLs(safeInt(entries, SvDataType.TELOMERIC_SGL))
                .build();
    }

    @NotNull
    private static List<CuppaPrediction> extractPredictions(@NotNull List<CuppaDataFile> files)
    {
        String bestCombinedType = determineBestCombinedDataType(files);
        if(bestCombinedType == null)
        {
            LOGGER.warn("Could not find a valid combined data type amongst cuppa entries");
            return Lists.newArrayList();
        }

        Map<String, CuppaDataFile> filesByType = files.stream()
                .filter(entry -> entry.Result.equals(ResultType.CLASSIFIER))
                .collect(Collectors.toMap(file -> file.DataType, file -> file));

        Map<ClassifierType, Map<String, Double>> predictionsByClassifier = Stream.of(ClassifierType.SNV_96_PAIRWISE,
                        ClassifierType.GENOMIC_POSITION_COHORT,
                        ClassifierType.FEATURE,
                        ClassifierType.ALT_SJ_COHORT,
                        ClassifierType.EXPRESSION_PAIRWISE)
                .collect(Collectors.toMap(classifier -> classifier, classifier -> predictionsForClassifier(filesByType, classifier)));

        return filesByType.get(bestCombinedType).CancerTypeValues.entrySet().stream().map(cancerPrediction ->
        {
            String cancerType = cancerPrediction.getKey();
            return ImmutableCuppaPrediction.builder()
                    .cancerType(cancerType)
                    .likelihood(cancerPrediction.getValue())
                    .snvPairwiseClassifier(predictionsByClassifier.get(ClassifierType.SNV_96_PAIRWISE).get(cancerType))
                    .genomicPositionClassifier(predictionsByClassifier.get(ClassifierType.GENOMIC_POSITION_COHORT).get(cancerType))
                    .featureClassifier(predictionsByClassifier.get(ClassifierType.FEATURE).get(cancerType))
                    .altSjCohortClassifier(predictionsByClassifier.get(ClassifierType.ALT_SJ_COHORT).get(cancerType))
                    .expressionPairwiseClassifier(predictionsByClassifier.get(ClassifierType.EXPRESSION_PAIRWISE).get(cancerType))
                    .build();
        }).sorted(new CuppaPredictionComparator()).collect(Collectors.toList());
    }

    @NotNull
    private static Map<String, Double> predictionsForClassifier(@NotNull Map<String, CuppaDataFile> filesByType,
            @NotNull ClassifierType classifierType)
    {
        CuppaDataFile cuppaDataFile = filesByType.get(classifierType.toString());
        return cuppaDataFile == null ? Collections.emptyMap() : cuppaDataFile.CancerTypeValues;
    }

    @Nullable
    private static String determineBestCombinedDataType(@NotNull List<CuppaDataFile> entries)
    {
        boolean hasDnaCombinedType = false;
        boolean hasRnaCombinedType = false;
        boolean hasOverallCombinedType = false;

        for(CuppaDataFile entry : entries)
        {
            if(entry.Category == CategoryType.COMBINED)
            {
                switch(entry.DataType)
                {
                    case DataTypes.DATA_TYPE_COMBINED:
                        hasOverallCombinedType = true;
                        break;
                    case DataTypes.DATA_TYPE_DNA_COMBINED:
                        hasDnaCombinedType = true;
                        break;
                    case DataTypes.DATA_TYPE_RNA_COMBINED:
                        hasRnaCombinedType = true;
                        break;
                    default:
                        LOGGER.warn("Unrecognized combined data type: {}", entry.DataType);
                        break;
                }
            }
        }

        if(hasOverallCombinedType)
        {
            return DataTypes.DATA_TYPE_COMBINED;
        }
        else if(hasDnaCombinedType)
        {
            return DataTypes.DATA_TYPE_DNA_COMBINED;
        }
        else if(hasRnaCombinedType)
        {
            return DataTypes.DATA_TYPE_RNA_COMBINED;
        }

        return null;
    }

    private static int safeInt(@NotNull List<CuppaDataFile> entries, @NotNull SvDataType svDataType)
    {
        CuppaDataFile entry = findSvEntry(entries, svDataType);
        if(entry != null)
        {
            return (int) Math.round(Double.parseDouble(entry.Value));
        }
        else
        {
            // -1 is a magic value that can never exist in reality.
            return -1;
        }
    }

    @Nullable
    private static CuppaDataFile findSvEntry(@NotNull List<CuppaDataFile> entries, @NotNull SvDataType dataType)
    {
        for(CuppaDataFile entry : entries)
        {
            if(entry.Category == CategoryType.SV && entry.DataType.equals(dataType.toString()))
            {
                return entry;
            }
        }

        LOGGER.warn("Could not find entry with data type '{}'", dataType);
        return null;
    }
}
