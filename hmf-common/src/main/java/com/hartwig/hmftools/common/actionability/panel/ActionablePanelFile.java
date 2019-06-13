package com.hartwig.hmftools.common.actionability.panel;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class ActionablePanelFile {

    private static final String DELIMITER = "\t";

    public static void write(@NotNull final String filename, @NotNull final List<ActionablePanel> panel) throws IOException {
        Files.write(new File(filename).toPath(), toLines(panel));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<ActionablePanel> panel) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        panel.stream().map(ActionablePanelFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "")
                .add("gene")
                .add("amplification")
                .add("deletion")
                .add("fusion")
                .add("variant")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final ActionablePanel panel) {
        return new StringJoiner(DELIMITER).add(String.valueOf(panel.gene()))
                .add(String.valueOf(panel.amplification()))
                .add(String.valueOf(panel.deletion()))
                .add(String.valueOf(panel.fusion()))
                .add(String.valueOf(panel.variant()))
                .toString();
    }

}
