package com.hartwig.hmftools.orange.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.exp;
import static net.sf.dynamicreports.report.builder.DynamicReports.margin;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import java.awt.Color;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.DynamicReports;
import net.sf.dynamicreports.report.builder.component.HorizontalListBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.constant.Evaluation;

public final class OrangeTemplateFactory
{
    public static VerticalListBuilder header(final String sampleId)
    {
        java.net.URL imageUrl = OrangeTemplateFactory.class.getResource("/orange_circos.png");
        if(imageUrl == null)
        {
            throw new IllegalStateException("Resource /orange_circos.png not found");
        }

        VerticalListBuilder imageContainer = cmp.verticalList(
                cmp.verticalGap(10),
                cmp.image(imageUrl.toString()).setFixedHeight(60)
        );

        VerticalListBuilder orangeBoxContainer = cmp.verticalList(
                cmp.text("SAMPLE").setStyle(OrangeFonts.SIDE_PANEL_LABEL_STYLE),
                cmp.text(sampleId).setStyle(OrangeFonts.SIDE_PANEL_VALUE_STYLE)
        ).setStyle(OrangeStyles.ORANGE_BOX_STYLE).setFixedDimension(170, 60);

        //        report.setPageMargin(margin().setTop(0).setRight(0).setLeft(20).setBottom(20));
        //
        //        return report
        //                .pageHeader(
        return cmp.verticalList(
                cmp.horizontalList(
                        cmp.horizontalGap(15),
                        imageContainer,
                        orangeHeaderTextComponent(),
                        orangeBoxContainer
                ),
                cmp.verticalGap(10)
        );
        //                );
        //                .pageFooter(DynamicReports.cmp.horizontalList(
        //                        DynamicReports.cmp.pageXofY().setStyle(OrangeFonts.PAGE_NUMBER_STYLE)
        //                ));
    }

    public static VerticalListBuilder footer()
    {
        return cmp.verticalList(
                cmp.verticalGap(30),
                cmp.pageXofY().setStyle(OrangeFonts.PAGE_NUMBER_STYLE)
        );
    }

    public static JasperReportBuilder applyHeaderLayout(final JasperReportBuilder report, final String sampleId)
    {
        java.net.URL imageUrl = OrangeTemplateFactory.class.getResource("/orange_circos.png");
        if(imageUrl == null)
        {
            throw new IllegalStateException("Resource /orange_circos.png not found");
        }

        VerticalListBuilder imageContainer = cmp.verticalList(
                cmp.verticalGap(10),
                cmp.image(imageUrl.toString()).setFixedHeight(60)
        );

        VerticalListBuilder orangeBoxContainer = cmp.verticalList(
                cmp.text("SAMPLE").setStyle(OrangeFonts.SIDE_PANEL_LABEL_STYLE),
                cmp.text(sampleId).setStyle(OrangeFonts.SIDE_PANEL_VALUE_STYLE)
        ).setStyle(OrangeStyles.ORANGE_BOX_STYLE).setFixedDimension(170, 60);

        report.setPageMargin(margin().setTop(0).setRight(0).setLeft(20).setBottom(20));

        return report
                .pageHeader(
                        cmp.horizontalList(
                                cmp.horizontalGap(15),
                                imageContainer,
                                orangeHeaderTextComponent(),
                                orangeBoxContainer
                        ),
                        cmp.verticalGap(10)
                );
        //                .pageFooter(DynamicReports.cmp.horizontalList(
        //                        DynamicReports.cmp.pageXofY().setStyle(OrangeFonts.PAGE_NUMBER_STYLE)
        //                ));
    }

    public static JasperReportBuilder applyFooter(final JasperReportBuilder report, final int totalPages)
    {
        return report.pageFooter(
                cmp.text(exp.jasperSyntax("\"\" + $V{PAGE_NUMBER} + \"/" + totalPages + "\"", String.class))
                        .setEvaluationTime(Evaluation.PAGE)
                        .setStyle(OrangeFonts.PAGE_NUMBER_STYLE)
                        .setFixedWidth(35)
        );
    }

    private static VerticalListBuilder orangeHeaderTextComponent()
    {
        Color[] letterColors = { OrangeColors.PALETTE_ORANGE_1, OrangeColors.PALETTE_ORANGE_2,
                OrangeColors.PALETTE_ORANGE_3, OrangeColors.PALETTE_ORANGE_4,
                OrangeColors.PALETTE_ORANGE_5, OrangeColors.PALETTE_ORANGE_6 };
        String[] letters = { "O", "R", "A", "N", "G", "E" };

        HorizontalListBuilder letterList = cmp.horizontalList();
        for(int i = 0; i < 6; i++)
        {
            letterList.add(cmp.text(letters[i]).setStyle(
                    stl.style().setFont(OrangeFonts.titleFont()).setPadding(0).setForegroundColor(letterColors[i])
            ));
        }

        HorizontalListBuilder fancyTextContainer = letterList.setFixedWidth(50);
        HorizontalListBuilder allTextContainer = cmp.horizontalList(
                cmp.filler(),
                fancyTextContainer,
                cmp.horizontalGap(1),
                cmp.text("Report").setStyle(
                        stl.style().setFont(OrangeFonts.titleFont()).setForegroundColor(OrangeColors.PALETTE_BLACK)
                ),
                cmp.filler()
        );

        return cmp.verticalList(cmp.verticalGap(30), allTextContainer);
    }

    private OrangeTemplateFactory()
    {
    }
}
