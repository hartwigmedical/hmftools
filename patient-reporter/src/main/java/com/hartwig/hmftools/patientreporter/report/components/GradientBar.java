package com.hartwig.hmftools.patientreporter.report.components;

import static net.sf.dynamicreports.report.builder.DynamicReports.cht;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;

import java.awt.Color;
import java.awt.GradientPaint;
import java.util.Optional;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import net.sf.dynamicreports.report.builder.column.TextColumnBuilder;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class },
             allParameters = true)
public abstract class GradientBar {
    @NotNull
    abstract Color startColor();

    @NotNull
    abstract Color endColor();

    @NotNull
    abstract String startText();

    @NotNull
    abstract String endText();

    @NotNull
    abstract Optional<Pair<Integer, String>> value();

    @NotNull
    abstract Optional<Pair<Integer, String>> marker();

    @NotNull
    public static GradientBar of(@NotNull final Color startColor, @NotNull final Color endColor, @NotNull final String startText,
            @NotNull final String endText, final int value) {
        return ImmutableGradientBar.builder()
                .startColor(startColor)
                .endColor(endColor)
                .startText(startText)
                .endText(endText)
                .value(Optional.of(ImmutablePair.of(value, Strings.EMPTY)))
                .build();
    }

    @NotNull
    public static GradientBar ofOnlyMarker(@NotNull final Color startColor, @NotNull final Color endColor, @NotNull final String startText,
            @NotNull final String endText,final int markerPosition) {
        return ImmutableGradientBar.builder()
                .startColor(startColor)
                .endColor(endColor)
                .startText(startText)
                .endText(endText)
                .marker(Optional.of(ImmutablePair.of(markerPosition, Strings.EMPTY)))
                .build();
    }

    @NotNull
    public static GradientBar of(@NotNull final Color startColor, @NotNull final Color endColor, @NotNull final String startText,
            @NotNull final String endText, final int value, final int markerPosition) {
        return ImmutableGradientBar.builder()
                .startColor(startColor)
                .endColor(endColor)
                .startText(startText)
                .endText(endText)
                .value(Optional.of(ImmutablePair.of(value, Strings.EMPTY)))
                .marker(Optional.of(ImmutablePair.of(markerPosition, Strings.EMPTY)))
                .build();
    }

    @NotNull
    public static GradientBar of(@NotNull final Color startColor, @NotNull final Color endColor, @NotNull final String startText,
            @NotNull final String endText) {
        return ImmutableGradientBar.builder()
                .startColor(startColor)
                .endColor(endColor)
                .startText(startText)
                .endText(endText)
                .build();
    }

    @NotNull
    public ComponentBuilder<?, ?> build() {
        final TextColumnBuilder<String> itemColumn = col.column("item", "item", String.class);
        final TextColumnBuilder<Integer> valueColumn = col.column("value", "value", Integer.class);
        final DRDataSource dataSource = new DRDataSource("item", "value");
        dataSource.add("value", 100);
        final GradientPaint gradientPaint = new GradientPaint(0, 0, startColor(), 10, 0, endColor());

        final GradientBarCustomizer customizer =
                ImmutableGradientBarCustomizer.of(gradientPaint, startText(), endText(), value(), marker());
        return cht.barChart()
                .customizers(customizer)
                .setCategory(itemColumn)
                .series(cht.serie(valueColumn))
                .setHeight(39)
                .setDataSource(dataSource);
    }
}
