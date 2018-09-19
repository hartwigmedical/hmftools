package com.hartwig.hmftools.actionability.compare_with_SOC;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.jetbrains.annotations.NotNull;

public class AnalyzerSOC {
    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(AnalyzerSOC.class);

    @NotNull
    private final List<AnalyzerSOC> methods;

    public AnalyzerSOC(@NotNull final List<AnalyzerSOC> methods) {
        this.methods = methods;
    }


    @NotNull
    public static AnalyzerSOC loadFilesMethods(String fileMethods) throws IOException {
        final List<AnalyzerSOC> methods = new ArrayList<>();

        return new AnalyzerSOC(methods);
    }
}
