package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
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
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.DoubleField;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.IntField;
import com.hartwig.hmftools.compar.common.field.LongField;
import com.hartwig.hmftools.compar.common.field.StringField;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public record IsofoxSummaryComparer(ComparConfig mConfig) implements ItemComparer
{
    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_SUMMARY;
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
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return List.of(
                new StringField(FLD_QC_STATUS, i -> IsofoxSummaryData.qcStatus(((IsofoxSummaryData) i).RnaStatistics().qcStatus()), true),
                new LongField(FLD_TOTAL_FRAGS, i -> ((IsofoxSummaryData) i).RnaStatistics().totalFragments(), true, 10., 0.01),
                new LongField(FLD_DUPLICATE_FRAGS, i -> ((IsofoxSummaryData) i).RnaStatistics().duplicateFragments(), true, 10., 0.01),
                new DoubleField(FLD_SPLICED_FRAG_PERC, i -> ((IsofoxSummaryData) i).RnaStatistics().splicedFragmentPerc(), true, 0.01, 0.05, "%.2f"),
                new DoubleField(FLD_UNSPLICED_FRAG_PERC, i -> ((IsofoxSummaryData) i).RnaStatistics().unsplicedFragmentPerc(), true, 0.01, 0.05, "%.2f"),
                new DoubleField(FLD_ALT_FRAG_PERC, i -> ((IsofoxSummaryData) i).RnaStatistics().altFragmentPerc(), true, 0.01, 0.05, "%.2f"),
                new DoubleField(FLD_CHIMERIC_FRAG_PERC, i -> ((IsofoxSummaryData) i).RnaStatistics().chimericFragmentPerc(), true, 0.01, 0.05, "%.2f"),
                new IntField(FLD_SPLICED_GENE_COUNT, i -> ((IsofoxSummaryData) i).RnaStatistics().splicedGeneCount(), true, 10., 0.01),
                new IntField(FLD_READ_LENGTH, i -> ((IsofoxSummaryData) i).RnaStatistics().readLength(), true, null, null),
                new DoubleField(FLD_FRAG_LENGTH_5TH, i -> ((IsofoxSummaryData) i).RnaStatistics().fragmentLength5thPercent(), true, null, 0.05, "%.1f"),
                new DoubleField(FLD_FRAG_LENGTH_50TH, i -> ((IsofoxSummaryData) i).RnaStatistics().fragmentLength50thPercent(), true, null, 0.05, "%.1f"),
                new DoubleField(FLD_FRAG_LENGTH_95TH, i -> ((IsofoxSummaryData) i).RnaStatistics().fragmentLength95thPercent(), true, null, 0.05, "%.1f"),
                new DoubleField(FLD_ENRICHED_GENE_PERC, i -> ((IsofoxSummaryData) i).RnaStatistics().enrichedGenePercent(), true, 0.01, null, "%.2f"),
                new DoubleField(FLD_MEDIAN_GC_RATIO, i -> ((IsofoxSummaryData) i).RnaStatistics().medianGCRatio(), true, 0.01, null, "%.2f"),
                new DoubleField(FLD_FORWARD_STRAND_PERC, i -> ((IsofoxSummaryData) i).RnaStatistics().forwardStrandPercent(), true, 0.01, null, "%.2f")
        );
    }

    @Override
    public List<String> displayFieldNames()
    {
        return List.of(
                FLD_QC_STATUS, FLD_TOTAL_FRAGS, FLD_DUPLICATE_FRAGS, FLD_SPLICED_FRAG_PERC, FLD_UNSPLICED_FRAG_PERC, FLD_ALT_FRAG_PERC,
                FLD_CHIMERIC_FRAG_PERC, FLD_READ_LENGTH, FLD_FRAG_LENGTH_5TH, FLD_FRAG_LENGTH_50TH, FLD_FRAG_LENGTH_95TH,
                FLD_ENRICHED_GENE_PERC
        );
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
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
        String filename = RnaStatisticFile.generateFilename(fileSources.Isofox, sampleId);
        String oldFilename = filename.replace(TSV_EXTENSION, CSV_EXTENSION);

        if(!Files.exists(Paths.get(filename)) && Files.exists(Paths.get(oldFilename)))
        {
            return oldFilename;
        }
        else
        {
            return filename;
        }
    }
}
