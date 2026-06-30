package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentileField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatTpmField;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.type;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.primitives.Doubles;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Expressions;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.jasperreports.engine.JREmptyDataSource;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

public final class ExpressionTable
{
    public static JasperReportBuilder build(final String title, final List<GeneExpression> expressions,
            final boolean sortAscending, final ReportResources reportResources)
    {
        if(expressions.isEmpty())
        {
            return report()
                    .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE),
                            cmp.text(ReportResources.NONE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE))
                    .setDataSource(new JREmptyDataSource());
        }

        List<Map<String, ?>> rows = new ArrayList<>();
        for(GeneExpression expression : sort(expressions, sortAscending))
        {
            String percentile;
            String foldChange;
            if(expression.percentileCancer() != null && expression.medianTpmCancer() != null)
            {
                percentile = formatPercentileField(expression.percentileCancer());
                foldChange = Expressions.formatFoldChangeCancer(expression);
            }
            else
            {
                percentile = formatPercentileField(expression.percentileCohort());
                foldChange = Expressions.formatFoldChange(expression);
            }

            Map<String, Object> row = new HashMap<>();
            row.put("gene", expression.gene());
            row.put("tpm", formatTpmField(expression.tpm()));
            row.put("percentile", percentile);
            row.put("foldchange", foldChange);
            row.put("empty", "");
            rows.add(row);
        }

        return report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE_WITH_GAP))
                .columns(
                        col.column("GENE", "gene", type.stringType()).setWidth(10),
                        col.column("TPM", "tpm", type.stringType()).setWidth(10),
                        col.column("PERCENTILE", "percentile", type.stringType()).setWidth(10),
                        col.column("FOLD CHANGE", "foldchange", type.stringType()).setWidth(10),
                        col.column("", "empty", type.stringType()).setWidth(30)
                )
                .setDataSource(new JRMapCollectionDataSource(rows));
    }

    private static List<GeneExpression> sort(final List<GeneExpression> expressions, final boolean sortAscending)
    {
        return expressions.stream().sorted((e1, e2) ->
        {
            if(sortAscending)
            {
                return Doubles.compare(e1.percentileCohort(), e2.percentileCohort());
            }
            else
            {
                return Doubles.compare(e2.percentileCohort(), e1.percentileCohort());
            }
        }).collect(Collectors.toList());
    }
}
