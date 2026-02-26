package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_ALT_FRAG_PERC;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_CHIMERIC_FRAG_PERC;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_DUPLICATE_FRAGS;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_ENRICHED_GENE_PERC;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_FORWARD_STRAND_PERC;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_FRAG_LENGTH_50TH;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_FRAG_LENGTH_5TH;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_FRAG_LENGTH_95TH;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_MEDIAN_GC_RATIO;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_QC_STATUS;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_READ_LENGTH;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_SPLICED_FRAG_PERC;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_SPLICED_GENE_COUNT;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_TOTAL_FRAGS;
import static com.hartwig.hmftools.compar.isofox.IsofoxSummaryData.FLD_UNSPLICED_FRAG_PERC;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.RnaStatisticFile;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public record IsofoxSummaryComparer(ComparConfig mConfig) implements ItemComparer
{
    @Override
    public CategoryType category()
    {
        return CategoryType.ISOFOX_SUMMARY;
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
        thresholds.addFieldThreshold(FLD_TOTAL_FRAGS, 10, 0.01);
        thresholds.addFieldThreshold(FLD_DUPLICATE_FRAGS, 10, 0.01);
        thresholds.addFieldThreshold(FLD_SPLICED_FRAG_PERC, 0.01, 0.05);
        thresholds.addFieldThreshold(FLD_UNSPLICED_FRAG_PERC, 0.01, 0.05);
        thresholds.addFieldThreshold(FLD_ALT_FRAG_PERC, 0.01, 0.05);
        thresholds.addFieldThreshold(FLD_CHIMERIC_FRAG_PERC, 0.01, 0.05);
        thresholds.addFieldThreshold(FLD_SPLICED_GENE_COUNT, 10, 0.01);
        thresholds.addFieldThreshold(FLD_FRAG_LENGTH_5TH, -1, 0.05);
        thresholds.addFieldThreshold(FLD_FRAG_LENGTH_50TH, -1, 0.05);
        thresholds.addFieldThreshold(FLD_FRAG_LENGTH_95TH, -1, 0.05);
        thresholds.addFieldThreshold(FLD_ENRICHED_GENE_PERC, 0.01, -1);
        thresholds.addFieldThreshold(FLD_MEDIAN_GC_RATIO, 0.01, -1);
        thresholds.addFieldThreshold(FLD_FORWARD_STRAND_PERC, 0.01, -1);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return List.of(
                FLD_QC_STATUS, FLD_TOTAL_FRAGS, FLD_DUPLICATE_FRAGS, FLD_SPLICED_FRAG_PERC, FLD_UNSPLICED_FRAG_PERC, FLD_ALT_FRAG_PERC,
                FLD_CHIMERIC_FRAG_PERC, FLD_READ_LENGTH, FLD_FRAG_LENGTH_5TH, FLD_FRAG_LENGTH_50TH, FLD_FRAG_LENGTH_95TH,
                FLD_ENRICHED_GENE_PERC
        );
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
            List<String> lines = Files.readAllLines(Paths.get(determineFileName(sampleId, fileSources)));
            RnaStatistics rnaStatistics = RnaStatisticFile.fromLines(lines);
            comparableItems.add(new IsofoxSummaryData(rnaStatistics));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Isofox Summary data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }

    private static String determineFileName(final String sampleId, final FileSources fileSources)
    {
        String current_file_name = RnaStatisticFile.generateFilename(fileSources.Isofox, sampleId);
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
