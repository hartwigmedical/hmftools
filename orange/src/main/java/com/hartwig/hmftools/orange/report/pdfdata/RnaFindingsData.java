package com.hartwig.hmftools.orange.report.pdfdata;

import java.util.List;

import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.datamodel.isofox.RnaStatistics;

public class RnaFindingsData
{
    public final boolean hasRnaFail;
    public final RnaStatistics summary;
    public final List<GeneExpression> highExpressionGenes;
    public final List<RnaFusion> fusions;
    public final List<NovelSpliceJunction> novelSpliceJunctions;

    public RnaFindingsData(
            final boolean hasRnaFail,
            final RnaStatistics summary,
            final List<GeneExpression> highExpressionGenes,
            final List<RnaFusion> fusions,
            final List<NovelSpliceJunction> novelSpliceJunctions)
    {
        this.hasRnaFail = hasRnaFail;
        this.summary = summary;
        this.highExpressionGenes = highExpressionGenes;
        this.fusions = fusions;
        this.novelSpliceJunctions = novelSpliceJunctions;
    }
}
