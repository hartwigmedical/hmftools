package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_ADJ_TPM;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.TranscriptExpressionFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public record IsofoxTranscriptDataComparer(ComparConfig mConfig) implements ItemComparer
{
    @Override
    public CategoryType category()
    {
        return CategoryType.ISOFOX_TRANSCRIPT_DATA;
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
        thresholds.addFieldThreshold(FLD_ADJ_TPM, -1, 0.05);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return List.of(FLD_GENE_NAME, FLD_ADJ_TPM);
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
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            String filename = determineFileName(sampleId, fileSources);
            TranscriptExpressionFile.read(filename).stream().map(IsofoxTranscriptData::new).forEach(comparableItems::add);
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Isofox Transcript data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }

    private static String determineFileName(final String sampleId, final FileSources fileSources)
    {
        String current_file_name = TranscriptExpressionFile.generateFilename(fileSources.Isofox, sampleId);
        String old_file_name = current_file_name.replace(".tsv", ".csv");

        if(!Files.exists(Paths.get(current_file_name)) && Files.exists(Paths.get(old_file_name)))
        {
            return old_file_name;
        }
        else
        {
            return current_file_name;
        }
    }
}
