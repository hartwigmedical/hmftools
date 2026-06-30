package com.hartwig.hmftools.orange.report.tables;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.type;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleChrArmCopyNumber;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.expression.JasperExpression;
import net.sf.jasperreports.engine.JREmptyDataSource;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

public final class ChrArmCopyNumberTable
{
    public static JasperReportBuilder build(
            final String title, final List<PurpleChrArmCopyNumber> chrArmCopyNumbers, final ReportResources reportResources)
    {
        if(chrArmCopyNumbers.isEmpty())
        {
            return report()
                    .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE),
                            cmp.text(ReportResources.NONE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE))
                    .setDataSource(new JREmptyDataSource());
        }

        List<Map<String, ?>> rows = new ArrayList<>();
        for(PurpleChrArmCopyNumber arm : sort(chrArmCopyNumbers))
        {
            Map<String, Object> row = new LinkedHashMap<>();
            row.put("chromosome", arm.chromosome());
            row.put("arm", arm.arm());
            row.put("type", arm.type());
            row.put("cn", formatSingleDigitDecimal(arm.copyNumber()));
            row.put("relcn", formatSingleDigitDecimal(arm.relativeCopyNumber()));
            row.put("driver", arm.driverInterpretation().toString());
            row.put("empty", "");
            rows.add(row);
        }

        return report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE_WITH_GAP))
                //                .pageHeader(
                //                        cmp.text("CONTINUED FROM THE PREVIOUS PAGE")
                //                                .setStyle(OrangeFonts.TABLE_CONTENT_STYLE)
                //                                .setPrintWhenExpression(new JasperExpression<>("$V{PAGE_NUMBER} > 1", Boolean.class))
                //                )
                .columns(
                        col.column("CHROMOSOME", "chromosome", type.stringType()).setWidth(10),
                        col.column("ARM", "arm", type.stringType()).setWidth(10),
                        col.column("TYPE", "type", type.stringType()).setWidth(10),
                        col.column("CN", "cn", type.stringType()).setWidth(10),
                        col.column("REL CN", "relcn", type.stringType()).setWidth(10),
                        col.column("DRIVER", "driver", type.stringType()).setWidth(10),
                        col.column("", "empty", type.stringType()).setWidth(30)
                )
                .pageFooter(cmp.text("THE TABLE CONTINUES ON THE NEXT PAGE").setStyle(OrangeFonts.TABLE_CONTENT_STYLE))
                .lastPageFooter(cmp.verticalGap(1))
                .setDataSource(new JRMapCollectionDataSource(rows));
    }

    private static List<PurpleChrArmCopyNumber> sort(final List<PurpleChrArmCopyNumber> arms)
    {
        return arms.stream().sorted((arm1, arm2) ->
        {
            DriverInterpretation d1 = arm1.driverInterpretation();
            DriverInterpretation d2 = arm2.driverInterpretation();
            if(d1 != d2)
            {
                return d1 == DriverInterpretation.HIGH ? -1 : 1;
            }
            int chr1 = HumanChromosome.chromosomeRank(arm1.chromosome());
            int chr2 = HumanChromosome.chromosomeRank(arm2.chromosome());
            if(chr1 != chr2)
            {
                return Integer.compare(chr1, chr2);
            }
            return arm1.arm().compareTo(arm2.arm());
        }).collect(Collectors.toList());
    }
}
