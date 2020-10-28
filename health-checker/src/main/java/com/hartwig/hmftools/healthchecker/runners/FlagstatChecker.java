package com.hartwig.hmftools.healthchecker.runners;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.healthchecker.result.ImmutableQCValue;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FlagstatChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(FlagstatChecker.class);

    @NotNull
    private final String refFlagstat;
    @Nullable
    private final String tumFlagstat;

    public FlagstatChecker(@NotNull final String refFlagstat, @Nullable final String tumFlagstat) {
        this.refFlagstat = refFlagstat;
        this.tumFlagstat = tumFlagstat;
    }

    public static double round(double value, int places) {
        if (places < 0) throw new IllegalArgumentException();
        BigDecimal bd = BigDecimal.valueOf(value);
        bd = bd.setScale(places, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }

    public static String divideTwoStrings(String string1, String string2) {
        Double double1 = Double.parseDouble(string1);
        Double double2 = Double.parseDouble(string2);
        Double proportion = round(double1 / double2, 6);
        return String.valueOf(proportion);
    }

    @Nullable
    public static String valueBySubstring(List<String> lines , String subString) {
        List<String> matchLines = new ArrayList<>();
        for(String line : lines) {
            if (line.contains(subString)) {
                matchLines.add(line);
            }
        }
        if (matchLines.size() == 1) {
            return matchLines.get(0).split(" ")[0];
        }
        return null;
    }

    @NotNull
    @Override
    public List<QCValue> run() throws IOException {

        List<QCValue> qcValues = Lists.newArrayList();

        List<String> refLines = Files.readAllLines(new File(refFlagstat).toPath());
        String refTotal = valueBySubstring( refLines, "total");
        String refMapped = valueBySubstring( refLines, "mapped (");
        String refDuplicates = valueBySubstring( refLines, "duplicates");
        if (refTotal == null || refMapped == null || refDuplicates == null){
            throw new IOException("Unable to parse flagstat file correctly");
        }
        String refMappingProportion = divideTwoStrings(refMapped, refTotal);
        String refDuplicateProportion = divideTwoStrings(refDuplicates, refTotal);
        qcValues.add(ImmutableQCValue.of(QCValueType.REF_PROPORTION_MAPPED, refMappingProportion));
        qcValues.add(ImmutableQCValue.of(QCValueType.REF_PROPORTION_DUPLICATE, refDuplicateProportion));

        if (tumFlagstat != null) {
            List<String> tumLines = Files.readAllLines(new File(tumFlagstat).toPath());
            String tumTotal = valueBySubstring( tumLines, "total");
            String tumMapped = valueBySubstring( tumLines, "mapped (");
            String tumDuplicates = valueBySubstring( tumLines, "duplicates");
            if (tumTotal == null || tumMapped == null || tumDuplicates == null){
                throw new IOException("Unable to parse flagstat file correctly");
            }
            String tumMappingProportion = divideTwoStrings(tumMapped, tumTotal);
            String tumDuplicateProportion = divideTwoStrings(tumDuplicates, tumTotal);
            qcValues.add(ImmutableQCValue.of(QCValueType.TUM_PROPORTION_MAPPED, tumMappingProportion));
            qcValues.add(ImmutableQCValue.of(QCValueType.TUM_PROPORTION_DUPLICATE, tumDuplicateProportion));
        }
        return qcValues;
    }
}
