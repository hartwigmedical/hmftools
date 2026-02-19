package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_ADJ_TPM;
import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_SPLICED_FRAGS;
import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_UNSPLICED_FRAGS;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.GeneExpressionFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public record IsofoxGeneDataComparer(ComparConfig mConfig) implements ItemComparer
{
    @Override
    public CategoryType category()
    {
        return CategoryType.ISOFOX_GENE_DATA;
    }

    @Override
    public boolean hasReportable()
    {
        return false;
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_SPLICED_FRAGS, 10, 0.05);
        thresholds.addFieldThreshold(FLD_UNSPLICED_FRAGS, 10, 0.05);
        thresholds.addFieldThreshold(FLD_ADJ_TPM, -1, 0.05);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return List.of(FLD_SPLICED_FRAGS, FLD_UNSPLICED_FRAGS, FLD_ADJ_TPM);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        // Not currently supported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        List<GeneExpression> geneExpressions = GeneExpressionFile.read(GeneExpressionFile.generateFilename(fileSources.Isofox, sampleId));
        if(geneExpressions == null){
            CMP_LOGGER.warn("sample({}) failed to load Isofox Gene data", sampleId);
            return null;
        }

        return geneExpressions.stream().<ComparableItem>map(IsofoxGeneData::new).toList();
    }
}
