package com.hartwig.hmftools.patientreporter.report.components;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import com.hartwig.hmftools.patientreporter.report.Commons;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.style.PenBuilder;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class },
             allParameters = true)
public abstract class SliderSection {
    private static final PenBuilder BORDER = stl.pen1Point();
    private static final int EDGE_PADDING = 7;
    private static final int INNER_PADDING = 2;

    @NotNull
    abstract String title();

    @NotNull
    abstract String details();

    @NotNull
    abstract String description();

    @NotNull
    abstract GradientBar bar();

    @NotNull
    public ComponentBuilder<?, ?> build() {
        return cmp.verticalList(dataSection(), descriptionSection());
    }

    @NotNull
    private ComponentBuilder<?, ?> dataSection() {
        return cmp.horizontalList(titleAndDetailsSection(),
                bar().build().setStyle(stl.style().setTopPadding(5).setRightPadding(5).setTopBorder(BORDER).setRightBorder(BORDER)));
    }

    @NotNull
    private ComponentBuilder<?, ?> titleAndDetailsSection() {
        return cmp.verticalList(cmp.text(title())
                        .setStyle(Commons.fontStyle()
                                .setFontSize(10)
                                .setBold(true)
                                .setLeftPadding(EDGE_PADDING)
                                .setTopPadding(EDGE_PADDING)
                                .setBottomPadding(INNER_PADDING)
                                .setTopBorder(BORDER)
                                .setLeftBorder(BORDER)),
                cmp.text(details())
                        .setStyle(Commons.monospaceFontStyle()
                                .setFontSize(10)
                                .setLeftPadding(EDGE_PADDING)
                                .setTopPadding(INNER_PADDING)
                                .setBottomPadding(INNER_PADDING)
                                .setLeftBorder(BORDER)));
    }

    @NotNull
    private ComponentBuilder<?, ?> descriptionSection() {
        return cmp.text(description())
                .setStyle(Commons.fontStyle()
                        .setFontSize(8)
                        .setTopPadding(INNER_PADDING)
                        .setLeftPadding(EDGE_PADDING)
                        .setRightPadding(EDGE_PADDING)
                        .setBottomPadding(EDGE_PADDING)
                        .setBottomBorder(BORDER)
                        .setLeftBorder(BORDER)
                        .setRightBorder(BORDER));
    }
}
