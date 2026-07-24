package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_AMP_DEL;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_CHROMOSOME_BAND;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonChromosome;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.DoubleField;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.StringField;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class GermlineAmpDelComparer implements ItemComparer
{
    protected static final String FLD_GERMLINE_STATUS = "GermlineStatus";
    protected static final String FLD_TUMOR_STATUS = "TumorStatus";
    protected static final String FLD_GERMLINE_CN = "GermlineCopyNumber";
    protected static final String FLD_TUMOR_CN = "TumorCopyNumber";

    private final ComparConfig mConfig;

    public GermlineAmpDelComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType category() { return GERMLINE_AMP_DEL; }

    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return List.of(
                new StringField(FLD_REPORTED, i -> ((GermlineAmpDelData) i).AmpDelData.Reported.toString(), true),
                new StringField(FLD_GERMLINE_STATUS, i -> ((GermlineAmpDelData) i).AmpDelData.NormalStatus.toString(),
                        true),
                new StringField(FLD_TUMOR_STATUS, i -> ((GermlineAmpDelData) i).AmpDelData.TumorStatus.toString(),
                        true),
                new DoubleField(FLD_GERMLINE_CN, i -> ((GermlineAmpDelData) i).AmpDelData.GermlineCopyNumber,
                        true, 0.2, 0.1, "%.2f"),
                new DoubleField(FLD_TUMOR_CN, i -> ((GermlineAmpDelData) i).AmpDelData.TumorCopyNumber,
                        true, 0.2, 0.1, "%.2f"),
                new StringField(FLD_CHROMOSOME, i -> ((GermlineAmpDelData) i).mComparisonChromosome,
                        true),
                new StringField(FLD_CHROMOSOME_BAND, i -> ((GermlineAmpDelData) i).AmpDelData.ChromosomeBand, true)
        );
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches, final FieldConfig fieldConfig)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches, fieldConfig);
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Lists.newArrayList(
                FLD_REPORTED, FLD_GERMLINE_STATUS, FLD_TUMOR_STATUS, FLD_GERMLINE_CN, FLD_TUMOR_CN);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
    {
        final List<GermlineAmpDel> germlineAmpDels = dbAccess.readGermlineCopyNumbers(sampleId);
        List<ComparableItem> items = Lists.newArrayList();
        germlineAmpDels.forEach(x -> items.add(createGermlineAmpDelData(x)));
        return items;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            String germlineAmpDelFile = GermlineAmpDel.generateFilename(fileSources.Purple, sampleId);

            if(!Files.exists(Paths.get(germlineAmpDelFile)))
            {
                // try pre v3.0 germline deletions
                germlineAmpDelFile = germlineAmpDelFile.replaceAll("purple.germline_amp_del.tsv", "purple.germline.deletion.tsv");
            }

            List<GermlineAmpDel> germlineAmpDels = GermlineAmpDel.read(germlineAmpDelFile);
            germlineAmpDels.forEach(x -> comparableItems.add(createGermlineAmpDelData(x)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to read germline amp-del data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }

    private GermlineAmpDelData createGermlineAmpDelData(final GermlineAmpDel deletion)
    {
        String comparisonChromosome = determineComparisonChromosome(deletion.Chromosome, mConfig.RequiresLiftover);
        return new GermlineAmpDelData(deletion, comparisonChromosome);
    }
}
